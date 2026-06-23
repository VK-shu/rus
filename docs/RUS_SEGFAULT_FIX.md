# rusSyn Core Dumped 修復紀錄

> **專案路徑：** `SliQSim/external/rus`
>
> **環境：** WSL (Ubuntu) on Windows
>
> **觸發指令：**
> ```bash
> ./rusSyn -T pi/8 -O out.qasm -F 20 -E 1e-3 -C g-count -D 2
> ```
>
> **修復後狀態：** `pi/8`、`pi/16` 均可正常產生 QASM（輸出於 `output/` 目錄）

---

## 1. 問題現象

執行 `rusSyn` 合成小角度 unitary 時，程式在 **Stage 3: Decomposition** 階段發生 segmentation fault (core dumped)。

典型輸出（修復前）：

```
Theta: 0.392699
Stage 1: PSLQ
input x:
(0.19509,0) (0.83147,0) (0.980785,0) (0.55557,0)
PSLQ: Iteration limit exceeded 200000
OmegaRing: -78+11*w +-21*w^2 +48*w^3
real epsilon = 1.72224e-08
Stage 2-1: Normalization
Stage 3-1: Decomposition
Segmentation fault (core dumped)
```

---

## 2. GDB 回溯分析

使用 GDB 確認崩潰點：

```
#0  resring<8>::operator<<=(int)
#1  ring_int<resring<8>>::gde2(resring<8>)
#2  ring_int_real<resring<8>>::gde() const
#3  exactDecomposer::decompose(...)
#4  RUS::run(...)
#5  main
```

**結論：** 崩潰發生在 `exactDecomposer::decompose()` 內部，呼叫 `ring_int<resring<8>>::gde2()` 時進入 `resring<8>::operator<<=`，觸發 SIGSEGV。

---

## 3. 根本原因（三點）

### 3.1 PSLQ 候選數量不足（主要原因）

| 項目 | 修復前 | 可運作的參考版本 (`final_project/rus`) |
|------|--------|------------------------------------------|
| PSLQ 迭代上限 | `10000 × effort`（`-F 20` → 200,000） | `run()` 內強制 `100,000,000` |
| 成功候選數 | 通常只有 1 個 `OmegaRing` | 可收集 10+ 個候選 |
| 後果 | 單一候選矩陣在 decompose 失敗即整體崩潰 | 多候選中總有可合成的 |

`-F 20` 代表「最多收集 20 個成功候選」，但內層 PSLQ 共用同一個 iteration 預算。預算太小時往往只成功 1 次，後續 Stage 3 沒有備選方案。

### 3.2 `resring` 運算子缺少 `return *this`（UB）

`core/resring.cpp` 中 `operator/=`、`operator<<=`、`operator>>=` 宣告回傳 reference，但實作未 `return *this`。

`gde2()` 迴圈內會執行 `a /= 2`，依賴運算子正確回傳值。此行為屬 C++ 未定義行為 (UB)，可能導致 stack 損壞，並在後續 `operator<<=` 處以 segfault 形式爆發。

### 3.3 PARI 重複初始化（潛在風險）

`rus/normalization.cpp` 的 `solve_y()` 每次呼叫都建立 `normSolver ns`，其建構子會執行 `pari_init_opts()`。`normsolver.h` 明確建議使用 `normSolver::instance()` 以避免多次初始化 PARI。

---

## 4. 所有保留的程式修改

> 以下為最終保留的改動。中間嘗試後又還原的項目（例如短暫改 iteration 預設為 `1000000 × effort`）不列入。

---

### 4.1 `rus/rus.cpp`

**① 建構子設定高精度**

```cpp
mpf_set_default_prec(2048);
mpfr_set_default_prec(2048);
```

- **原因：** 小角度合成需要足夠精度，與可運作參考版本一致。

**② `run()` 內 PSLQ 迭代上限（僅在未指定 `-I` 時套用）**

```cpp
if (!pslq_iters_from_cli)
    pslq_max_iter = DEFAULT_PSLQ_MAX_ITER;  // 100000000
```

- **原因：** 最關鍵修復。預設確保 PSLQ 有足夠迭代預算收集多個候選；若 CLI 有 `-I` 則尊重使用者設定。

**③ 移除 3 參數建構子**

- 只保留 4 參數建構子：`RUS(debug, epsilon, effort, pslq_iters, criterion)`
- **原因：** iteration 預設政策統一由 `apps/main_rus.cpp` 管理，避免兩處各自定義造成混淆。

---

### 4.2 `rus/rus.h`

- 移除 3 參數建構子宣告
- 加入 `#include <mpfr.h>`（配合高精度設定）

---

### 4.3 `apps/main_rus.cpp`

**`-I` / iteration 行為**

- 未指定 `-I`：傳入 `pslq_iters_from_cli = false`，由 `run()` 使用 `DEFAULT_PSLQ_MAX_ITER`（100,000,000）
- 指定 `-I <N>`：傳入 `pslq_iters_from_cli = true`，`run()` 保留 CLI 的 `N`，不再覆寫

```cpp
const bool pslq_iters_from_cli = vm.count("iterations") > 0;
RUS rus(debug_level, epsilon, effort, iterations, criterion, pslq_iters_from_cli);
```

- **原因：** 預設維持高迭代上限以確保穩定合成，同時允許實驗時用 `-I` 手動調整。

**輸出路徑（後續整理）**

- 相對路徑的 `-O` 預設寫入 `output/` 目錄，避免 QASM 與原始碼混在同一層。

### 4.4 `rus/normalization.cpp`

**修改前：**

```cpp
normSolver ns;
```

**修改後：**

```cpp
const normSolver &ns = normSolver::instance();
```

- **原因：** 避免重複 `pari_init`，與專案其他模組（`core/factorzs2.cpp`、`appr/toptzrot2.cpp` 等）用法一致。

---

### 4.5 `core/resring.cpp`

為以下運算子補上 `return *this;`：

- `operator/=`
- `operator<<=`
- `operator>>=`

**修改前（以 `operator<<=` 為例）：**

```cpp
resring<TMod> &resring<TMod>::operator <<=(int e)
{
    set( v << e );
}
```

**修改後：**

```cpp
resring<TMod> &resring<TMod>::operator <<=(int e)
{
    set( v << e );
    return *this;
}
```

- **原因：** 修復 C++ UB；`gde2()` 內 `a /= 2` 依賴正確回傳 reference。

---

### 4.6 `CMakeLists.txt`

- CMake 最低版本提升至 3.16
- 啟用 C++17（`CMAKE_CXX_STANDARD 17`）
- 改用 target-based 連結（`target_link_libraries`）
- 透過 `find_path` / `find_library` + pkg-config 尋找 GMP、MPFR、PARI
- 原始碼路徑更新為 `apps/`、`core/` 等子目錄

- **原因：** 在 WSL/Linux 環境正確找到並連結依賴，確保可成功編譯 `rusSyn`。

---

### 4.7 `README.md`

- 更新依賴清單（C++17、Boost、GMP、MPFR、PARI）
- 補充 WSL/Ubuntu 安裝與編譯步驟
- 補充 PARI 從原始碼安裝與 `~/.gprc` 設定說明

- **原因：** 降低環境設定失敗機率（PARI 缺失或 `gprc` 未設定是常見建置問題）。

---

## 5. 環境設定（非程式碼，但必要）

在 WSL 中若 PARI 執行期找不到設定檔，需執行：

```bash
cp /etc/gprc ~/.gprc
```

此步驟本身**不能**單獨解決 segfault，但缺少時可能導致 norm solver 異常。

---

## 6. effort 與 iteration 的關係

| 參數 | CLI | 作用 |
|------|-----|------|
| **effort** | `-F` | 外層迴圈：最多收集幾個 PSLQ 成功候選 |
| **iterations** | `-I` | 內層迴圈：PSLQ 總迭代預算（所有 effort 輪共用） |

**預設（未給 `-I`）：**

```
pslq_max_iter = DEFAULT_PSLQ_MAX_ITER  // 100,000,000，在 run() 內設定
```

**有給 `-I <N>`：**

```
pslq_max_iter = N  // 使用 CLI 值，run() 不覆寫
```

**PSLQ 雙層迴圈邏輯（`pslq.cpp`）：**

```
it_count = itm   // 只初始化一次，不會每輪 effort 重置

while (effort 輪數 > 0):
    while (未成功 and it_count > 0):
        pslq_iteration()
    if 未成功:
        break
    將候選加入 r_sets
    重置 status，繼續下一輪 effort
```

**目前實際行為：** 未指定 `-I` 時 `run()` 使用 1 億次上限；指定 `-I` 時由 CLI 控制。`-F` 仍控制最多收集幾個候選。

---

## 7. 驗證方式

### 編譯

```bash
cd external/rus
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

或使用 in-source build：

```bash
cmake -DCMAKE_BUILD_TYPE=Release .
make -j$(nproc)
```

### 執行測試

```bash
./rusSyn -T pi/8 -O out.qasm -F 20 -E 1e-3 -C g-count -D 2
./rusSyn -T pi/16 -O out_pi16.qasm -F 20 -E 1e-3 -C g-count -D 1
```

### 預期結果（修復後）

- 不再 segfault
- 可看到多輪 `Stage 2-N` / `Stage 3-N`（通常 10+ 個候選）
- 成功寫出 `output/out.qasm`
- 可能仍印出 `PSLQ: Iteration limit exceeded 100000000`，但不影響最終結果

---

## 8. 修改檔案總覽

| 檔案 | 修改類型 |
|------|----------|
| `rus/rus.cpp` | 高精度、PSLQ 迭代上限、建構子簡化 |
| `rus/rus.h` | 建構子宣告、mpfr include |
| `apps/main_rus.cpp` | iteration 預設公式、輸出至 `output/` |
| `rus/normalization.cpp` | `normSolver::instance()` |
| `core/resring.cpp` | 運算子 `return *this` |
| `CMakeLists.txt` | 現代化建置與依賴尋找、新路徑 |
| `README.md` | 文件更新 |

---

## 9. 參考對照

可運作版本位於：

```
school/QDA/final_project/rus
```

關鍵差異為 `rus.cpp` 中 `run()` 的 `pslq_max_iter = 100000000`，以及建構子內 2048 位精度設定。其餘演算法結構與 `SliQSim/external/rus` 相同。

---

## 10. 已知限制與後續可選改善

- `PSLQ: Iteration limit exceeded` 訊息仍可能出現，屬正常現象
- README 註記：高 effort / precision 下程式仍可能不穩定（原 upstream 限制）
