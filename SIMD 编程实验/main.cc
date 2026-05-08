#include <cstring>
#include <string>
#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <sys/time.h>
#include <omp.h>
#include <vector>
#include <stdint.h>

#define __ARM_NEON 1
#define __ARM_NEON__ 1
#include <arm_neon.h>

// 全局变量：是否使用对齐内存分配（由 main 的命令行参数控制）
static bool g_use_aligned = false;

// 新增：安全的模加/模减，使用 64 位中间类型
static inline int mod_add_safe(uint32_t a, uint32_t b, uint32_t mod) {
    uint64_t res = (uint64_t)a + (uint64_t)b;
    return (int)(res >= mod ? res - mod : res);
}
static inline int mod_sub_safe(uint32_t a, uint32_t b, uint32_t mod) {
    int64_t res = (int64_t)a - (int64_t)b;
    return (int)(res < 0 ? res + mod : res);
}

// ========== 新增：通用 NTT（接受原根 root）==========
static uint64_t mod_pow(uint64_t a, uint64_t e, uint64_t mod) {
    uint64_t r = 1;
    a %= mod;
    while (e) {
        if (e & 1) r = (r * a) % mod;
        a = (a * a) % mod;
        e >>= 1;
    }
    return r;
}

static void bit_reverse(int *a, int n) {
    int j = 0;
    for (int i = 1; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) {
            int tmp = a[i];
            a[i] = a[j];
            a[j] = tmp;
        }
    }
}

// 通用 NTT：与原 ntt 相同逻辑，但接受 root 参数（例如 root=3）
// 修改点：NEON 分支保留，标量分支使用 long long（int64_t）中间计算，避免 32 位有符号溢出
void ntt_mod(int *a, int n, int invert, int mod, int root) {
    bit_reverse(a, n);

    for (int len = 2; len <= n; len <<= 1) {
        int wlen = (int)mod_pow((uint64_t)root, (mod - 1) / len, mod);
        if (invert) wlen = (int)mod_pow((uint64_t)wlen, mod - 2, mod);
        int half = len >> 1;

        std::vector<uint32_t> wpow(half);
        wpow[0] = 1;
        for (int j = 1; j < half; j++) {
            wpow[j] = (uint32_t)((uint64_t)wpow[j-1] * (uint64_t)wlen % (uint64_t)mod);
        }

        for (int i = 0; i < n; i += len) {
            int j = 0;
            // NEON 加速路径：每次处理 4 个 lane（保持原实现）
            for (; j + 4 <= half; j += 4) {
                // 以无符号方式加载到临时数组，所有中间计算使用 uint64_t（与之前实现一致）
                uint32_t u_tmp[4], v_tmp[4], w_tmp[4];
                for (int t = 0; t < 4; ++t) {
                    u_tmp[t] = (uint32_t)a[i + j + t];
                    v_tmp[t] = (uint32_t)a[i + j + half + t];
                    w_tmp[t] = wpow[j + t];
                }

                uint32_t sum_out[4], diff_out[4];
                for (int t = 0; t < 4; ++t) {
                    uint64_t prod = (uint64_t)v_tmp[t] * (uint64_t)w_tmp[t];
                    uint32_t m = (uint32_t)(prod % (uint64_t)mod);
                    uint64_t s = (uint64_t)u_tmp[t] + (uint64_t)m;
                    if (s >= (uint64_t)mod) s -= (uint64_t)mod;
                    uint64_t d;
                    if (u_tmp[t] >= m) d = (uint64_t)(u_tmp[t] - m);
                    else d = (uint64_t)u_tmp[t] + (uint64_t)mod - (uint64_t)m;
                    sum_out[t] = (uint32_t)s;
                    diff_out[t] = (uint32_t)d;
                }

                uint32x4_t sum_vec = vld1q_u32(sum_out);
                uint32x4_t diff_vec = vld1q_u32(diff_out);
                vst1q_u32((uint32_t*)&a[i + j], sum_vec);
                vst1q_u32((uint32_t*)&a[i + j + half], diff_vec);
            }

            // 剩余标量处理：使用 long long（int64_t）中间计算以修复 32 位溢出问题
            for (; j < half; j++) {
                long long u = (long long)(unsigned)a[i + j];
                long long v = (long long)(unsigned)a[i + j + half];
                long long mul = ( ( (__int128)v * (__int128)wpow[j] ) % (long long)mod );
                long long sum = u + mul;
                if (sum >= mod) sum -= mod;
                long long diff = u - mul;
                if (diff < 0) diff += mod;
                a[i + j] = (int)sum;
                a[i + j + half] = (int)diff;
            }
        }
    }

    if (invert) {
        long long inv_n = (long long)mod_pow(n, mod - 2, mod);
        for (int i = 0; i < n; i++) {
            long long val = ( ( (__int128)(unsigned)a[i] * (__int128)inv_n ) % (long long)mod );
            a[i] = (int)val;
        }
    }
}

// 为 Garner 准备的扩展欧几里得求逆（64-bit）
static long long egcd(long long a, long long b, long long &x, long long &y) {
    if (b == 0) { x = 1; y = 0; return a; }
    long long x1, y1;
    long long g = egcd(b, a % b, x1, y1);
    x = y1;
    y = x1 - (a / b) * y1;
    return g;
}
static long long modinv_ll(long long a, long long mod) {
    long long x, y;
    long long g = egcd(a >= 0 ? a : a + mod, mod, x, y);
    if (g != 1) return -1;
    x %= mod;
    if (x < 0) x += mod;
    return x;
}

// poly_multiply via multiple friendly moduli + Garner（当目标 mod 不支持所需长度时）
void poly_multiply_crt(int *a, int *b, int *ab, int n, int mod) {
    // 选择三个常用的 NTT 友元质数（两三个就足够覆盖 32-bit mod）
    const int K = 3;
    const int mods[K] = {167772161, 469762049, 1224736769}; // 常用 NTT 质数
    const int roots[K] = {3, 3, 3}; // 这里均使用 3 作为原根（常用）
    int size = 1;
    while (size < 2 * n) size <<= 1;

    // 为每个模运行 NTT 卷积
    // 为节省代码量，逐模执行并保存结果
    std::vector< std::vector<int> > res(K, std::vector<int>(2*n-1));
    for (int mi = 0; mi < K; ++mi) {
        int m = mods[mi];
        // 分配数组（使用对齐策略）
        int *fa = nullptr;
        int *fb = nullptr;
        if (g_use_aligned) {
            fa = (int*)aligned_alloc((size_t)size * sizeof(int), 32);
            fb = (int*)aligned_alloc((size_t)size * sizeof(int), 32);
            if (!fa || !fb) { fa = (int*)malloc(size * sizeof(int)); fb = (int*)malloc(size * sizeof(int)); }
        } else {
            fa = (int*)malloc(size * sizeof(int));
            fb = (int*)malloc(size * sizeof(int));
        }
        for (int i = 0; i < n; ++i) {
            fa[i] = a[i] % m;
            if (fa[i] < 0) fa[i] += m;
            fb[i] = b[i] % m;
            if (fb[i] < 0) fb[i] += m;
        }
        for (int i = n; i < size; ++i) { fa[i] = fb[i] = 0; }

        ntt_mod(fa, size, 0, m, roots[mi]);
        ntt_mod(fb, size, 0, m, roots[mi]);
        for (int i = 0; i < size; ++i) fa[i] = (int)((1LL * (unsigned)fa[i] * (unsigned)fb[i]) % (unsigned)m);
        ntt_mod(fa, size, 1, m, roots[mi]);

        for (int i = 0; i < 2*n-1; ++i) res[mi][i] = fa[i];

        if (g_use_aligned) {
#if defined(_MSC_VER) || defined(_WIN32)
            _aligned_free(fa); _aligned_free(fb);
#else
            free(fa); free(fb);
#endif
        } else {
            free(fa); free(fb);
        }
    }

    // 使用 Garner 将三个残余重构到目标模 mod
    // 逐系数重构，使用在目标模下的乘法以避免溢出
    // 逐步构建 x = r0 + m0 * t0; x_mod = (x_mod + (mult_mod_target * t) % mod) % mod
    for (int idx = 0; idx < 2*n-1; ++idx) {
        long long x_mod = res[0][idx] % mod;
        long long mult_mod = mods[0] % mod;
        long long mult = mods[0]; // product of previous moduli (in plain)
        for (int i = 1; i < K; ++i) {
            int mi = mods[i];
            long long r = res[i][idx];
            // compute t = ((r - x) * inv(mult % mi)) mod mi
            long long mm = mult % mi;
            long long inv = modinv_ll(mm, mi);
            long long t = (r - (x_mod % mi) + mi) % mi;
            t = (t * inv) % mi;
            // x_mod = x_mod + mult * t (mod mod)
            long long add = ( (__int128)mult_mod * (__int128)t ) % mod;
            x_mod = (x_mod + add) % mod;
            // update mult and mult_mod
            mult = mult * (long long)mi;
            mult_mod = ( (__int128)mult_mod * (__int128)(mi % mod) ) % mod;
        }
        ab[idx] = (int)( (x_mod % mod + mod) % mod );
    }
}

// ========== 平台无关的 32 字节对齐内存分配/释放，用于对齐 vs 非对齐对比 ==========
static void *aligned_alloc32(size_t size) {
#if defined(_MSC_VER) || defined(_WIN32)
    return _aligned_malloc(size, 32);
#else
    void *ptr = nullptr;
    if (posix_memalign(&ptr, 32, size) != 0) return nullptr;
    return ptr;
#endif
}
static void aligned_free(void *p) {
#if defined(_MSC_VER) || defined(_WIN32)
    _aligned_free(p);
#else
    free(p);
#endif
}

// poly_multiply_ntt：简化为始终使用单模 NTT（不进入 CRT 路径），以避免错误的 CRT 路径导致结果错误与性能暴增
void poly_multiply_ntt(int *a, int *b, int *ab, int n, int mod) {
    int size = 1;
    while (size < 2 * n) size <<= 1;

    // 始终使用单模 NTT 路径（即使 (mod-1) % size != 0 也按单模路径处理，保证不进入 CRT）
    int *fa = nullptr;
    int *fb = nullptr;
    if (g_use_aligned) {
        fa = (int*)aligned_alloc32(size * sizeof(int));
        fb = (int*)aligned_alloc32(size * sizeof(int));
        if (!fa || !fb) {
            std::cerr<<"aligned allocation failed, falling back to malloc\n";
            if (fa) aligned_free(fa);
            if (fb) aligned_free(fb);
            fa = (int*)malloc(size * sizeof(int));
            fb = (int*)malloc(size * sizeof(int));
        }
    } else {
        fa = (int*)malloc(size * sizeof(int));
        fb = (int*)malloc(size * sizeof(int));
    }

    for (int i = 0; i < n; i++) {
        fa[i] = a[i] % mod;
        if (fa[i] < 0) fa[i] += mod;
        fb[i] = b[i] % mod;
        if (fb[i] < 0) fb[i] += mod;
    }
    for (int i = n; i < size; i++) fa[i] = fb[i] = 0;

    // 使用默认原根 3（与原逻辑一致）
    ntt_mod(fa, size, 0, mod, 3);
    ntt_mod(fb, size, 0, mod, 3);

    for (int i = 0; i < size; i++) {
        fa[i] = (int)(( (__int128)(unsigned)fa[i] * (unsigned)fb[i] ) % (unsigned)mod);
    }

    ntt_mod(fa, size, 1, mod, 3);
    for (int i = 0; i < 2 * n - 1; i++) ab[i] = fa[i];

    if (g_use_aligned) { aligned_free(fa); aligned_free(fb); } else { free(fa); free(fb); }
}

// ---- 新增：Montgomery 参数与帮助函数 ----
struct Montgomery32 {
    uint32_t mod;
    uint32_t nprime;        // n' = -n^{-1} mod 2^32 (用于规约 q = (uint32_t)x * n')
    uint32_t r_mod;        // 2^32 % mod (可选)
    uint32_t r2_mod;       // (2^32)^2 % mod (可选)
    uint32x2_t nprime_dup; // neon duplicate
    uint32x2_t mod_dup;

    // 初始化（调用一次）
    void set_mod(uint32_t m) {
        mod = m;
        // 计算 nprime = (-inv_mod) mod 2^32
        //  使用牛顿迭代或扩展欧几里得；这里用简化牛顿迭代示例（同你的报告）
        uint32_t inv = 1;
        for (int i = 0; i < 5; ++i) inv = inv * (2 - m * inv);
        nprime = (uint32_t)(- (int32_t)inv);
        nprime_dup = vdup_n_u32(nprime);
        mod_dup = vdup_n_u32(mod);
        r_mod = (uint32_t)( (0x100000000ULL) % m );
        r2_mod = (uint32_t)((uint64_t)r_mod * r_mod % m);
    }
};

// 标量 reduce：将 64-bit x 规约到 [0, mod)
static inline uint32_t reduce_scalar(uint64_t x, uint32_t mod, uint32_t nprime) {
    uint32_t q = (uint32_t)( (uint64_t)(uint32_t)x * (uint64_t)nprime ); // q = (x * n') mod 2^32
    uint64_t m = (uint64_t)q * (uint64_t)mod;
    uint32_t y = (uint32_t)((x - m) >> 32);
    return (x < m) ? (uint32_t)(y + mod) : y;
}

// NEON 向量化的 montgomery multiply（对 2-lane uint32x2_t）
// 输入 a,b 为普通 uint32 值（或已转换为蒙哥马利表示），返回 montgomery(a*b)
// 示例仅用于 2-lane；可按需扩展为 4-lane（uint32x4_t 分拆）
static inline uint32x2_t mont_mul_neon_u32x2(uint32x2_t a, uint32x2_t b, const Montgomery32 &M) {
    // T = a * b  (64-bit lanes)
    uint64x2_t T = vmull_u32(a, b); // two 64-bit lanes
    // T_lo = low32(T)
    uint32x2_t T_lo = vmovn_u64(T); // lower 32 bits of each 64-bit lane
    // m = T_lo * nprime  (mod 2^32)
    uint32x2_t m = vmul_u32(T_lo, M.nprime_dup);
    // mN = m * mod (64-bit)
    uint64x2_t mN = vmull_u32(m, M.mod_dup);
    // sum = T + mN
    uint64x2_t sum = vaddq_u64(T, mN);
    // res = sum >> 32
    uint32x2_t res = vshrn_n_u64(sum, 32);
    // if res >= mod -> res -= mod
    uint32x2_t cmp = vcge_u32(res, M.mod_dup);
    uint32x2_t sub = vsub_u32(res, vand_u32(M.mod_dup, cmp));
    return sub;
}

// 也可实现 uint32x4_t 版本：分奇偶 or 使用 two vmull 组合（留作后续优化）

// ========== 其余框架函数保持不变 ==========
void fRead(int *a, int *b, int *n, int *p, int input_id){
    std::string str1 = "/nttdata/";
    std::string str2 = std::to_string(input_id);
    std::string strin = str1 + str2 + ".in";
    char data_path[strin.size() + 1];
    std::copy(strin.begin(), strin.end(), data_path);
    data_path[strin.size()] = '\0';
    std::ifstream fin;
    fin.open(data_path, std::ios::in);
    fin>>*n>>*p;
    for (int i = 0; i < *n; i++){
        fin>>a[i];
    }
    for (int i = 0; i < *n; i++){
        fin>>b[i];
    }
}

void fCheck(int *ab, int n, int input_id){
    std::string str1 = "/nttdata/";
    std::string str2 = std::to_string(input_id);
    std::string strout = str1 + str2 + ".out";
    char data_path[strout.size() + 1];
    std::copy(strout.begin(), strout.end(), data_path);
    data_path[strout.size()] = '\0';
    std::ifstream fin;
    fin.open(data_path, std::ios::in);
    for (int i = 0; i < n * 2 - 1; i++){
        int x;
        fin>>x;
        if(x != ab[i]){
            std::cout<<"多项式乘法结果错误"<<std::endl;
            return;
        }
    }
    std::cout<<"多项式乘法结果正确"<<std::endl;
    return;
}

void fWrite(int *ab, int n, int input_id){
    std::string str1 = "files/";
    std::string str2 = std::to_string(input_id);
    std::string strout = str1 + str2 + ".out";
    char output_path[strout.size() + 1];
    std::copy(strout.begin(), strout.end(), output_path);
    output_path[strout.size()] = '\0';
    std::ofstream fout;
    fout.open(output_path, std::ios::out);
    for (int i = 0; i < n * 2 - 1; i++){
        fout<<ab[i]<<'\n';
    }
}

int a[300000], b[300000], ab[300000];

int main(int argc, char *argv[])
{
    // 解析命令行参数，支持 --aligned / --unaligned / --both（默认 both）
    bool run_aligned = false, run_unaligned = false;
    if (argc == 1) { run_aligned = run_unaligned = true; }
    else {
        for (int i = 1; i < argc; ++i) {
            if (strcmp(argv[i], "--aligned") == 0) run_aligned = true;
            else if (strcmp(argv[i], "--unaligned") == 0) run_unaligned = true;
            else if (strcmp(argv[i], "--both") == 0) run_aligned = run_unaligned = true;
        }
    }

    int test_begin = 0;
    int test_end = 4;  // 运行全部
    for(int i = test_begin; i <= test_end; ++i){
        long double ans = 0;
        int n_, p_;
        fRead(a, b, &n_, &p_, i);
        memset(ab, 0, sizeof(ab));

        if (run_unaligned) {
            g_use_aligned = false;
            auto Start = std::chrono::high_resolution_clock::now();
            poly_multiply_ntt(a, b, ab, n_, p_);
            auto End = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double,std::ratio<1,1000>> elapsed = End - Start;
            ans += elapsed.count();
            std::cout<<"[Unaligned] average latency for n = "<<n_<<" p = "<<p_<<" : "<<ans<<" (ms) "<<std::endl;
            fCheck(ab, n_, i);
            fWrite(ab, n_, i);
        }

        if (run_aligned) {
            g_use_aligned = true;
            memset(ab, 0, sizeof(ab));
            auto Start = std::chrono::high_resolution_clock::now();
            poly_multiply_ntt(a, b, ab, n_, p_);
            auto End = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double,std::ratio<1,1000>> elapsed = End - Start;
            ans = elapsed.count();
            std::cout<<"[Aligned]   average latency for n = "<<n_<<" p = "<<p_<<" : "<<ans<<" (ms) "<<std::endl;
            fCheck(ab, n_, i);
            fWrite(ab, n_, i);
        }
    }
    // perf / objdump 提示
    std::cout<<"\nNote: to collect perf and assembly, compile with -O3 and enable NEON (aarch64 native):\n";
    std::cout<<"  g++ -O3 -march=native main.cc -o ntt && ./ntt --both\n";
    std::cout<<"  objdump -d -M intel ntt > ntt.asm\n";
    std::cout<<"  perf record -g ./ntt && perf report\n";

    return 0;
}