# 3DAlphaWrap

Standalone C++17 implementation of *Alpha Wrapping with an Offset* (Portaneri et al., ACM TOG 2022) — no external libraries, no CGAL.

对论文 *Alpha Wrapping with an Offset*（Portaneri et al., ACM TOG 2022）的一个**零依赖** C++17 实现：
输入 OBJ 三角面片汤（可以有点、洞、自交、翻面），输出**严格包住输入、水密、可定向**的三角网格 OBJ。

不依赖 CGAL / Eigen / Boost / OpenMesh，只用 C++ 标准库。单文件、单线程，约 1200 行。

## 1. 编译

MinGW / g++：

```bat
g++ -O2 -std=c++17 -Wall -Wextra -o alpha_wrap.exe src\alpha_wrap.cpp
```

CMake：

```bat
cmake -S . -B build -G "Visual Studio 17 2022" -A x64
cmake --build build --config Release
```

MSVC（在 “x64 Native Tools Command Prompt” 中）：

```bat
cl /nologo /O2 /std:c++17 /EHsc src\alpha_wrap.cpp /Fe:alpha_wrap.exe
```

## 2. 用法

```bat
alpha_wrap input.obj output.obj [--alpha A] [--offset D] [--verbose] [--debug]
```

| 参数 | 含义 |
|---|---|
| `--alpha A` | 空球最小半径（绝对长度）。默认 `0.02 * bbox 对角线` |
| `--offset D` | 偏移距离 δ（绝对长度）。默认 `0.002 * bbox 对角线` |
| `--verbose` | 额外打印输出顶点到输入的距离统计 |
| `--debug` | 打印逐次插入/腔体/门诊断（开发用） |

示例：

```bat
alpha_wrap test\sphere.obj test\sphere_wrap.obj --alpha 0.04 --offset 0.004
alpha_wrap test\box.obj    test\box_wrap.obj    --alpha 0.1  --offset 0.01
```

## 3. 算法对应关系

代码按论文 Algorithm 1 组织：

| 论文 | 代码 |
|---|---|
| 初始化：偏移面包围盒 + 8 角点 Delaunay | `main()` 中 `BOXT` / `makeTet` / 链接循环 |
| 无限单元（凸包外） | 用“哨兵邻居 `nb = -1`”表示：凸包面片的另一侧即外部 |
| gate 与 Φ(f) | `gatePhi()`：沿对偶边取最小 Delaunay 球半径 |
| 优先队列（按 Φ 降序） | `std::priority_queue<Gate>` |
| R1：对偶 Voronoi 边 × 偏移面 | `segmentCrossesOffset()`（二分） |
| R2：外接球心 → 输入投影 → 拉回偏移面 | `offsetPointTowardCenter()`（沿 [k, proj] 二分） |
| 雕刻 | `state = 1`，并把内侧邻接面片入队 |
| 输出 inside/outside 分界面 | 末尾遍历存活单元，输出 `FACE` 定向的面片 |

关键实现点：

- **Delaunay**：增量 Bowyer–Watson（`insertPoint`）。先定位包含点的四面体，再从它出发 BFS 出“外接球包含新点”的冲突区；若边界不是星形（数值退化），抖动重试。
- **谓词**：`orient3d` / `insphere` 用 `long double` + 相对误差过滤；MSVC 下 `long double == double`，用 `1e-12` 相对阈值。不是精确谓词，见“限制”。
- **距离 oracle**：三角面片 BVH（`buildBVH`），支持最近点/最近距离（`closestToInput` / `distToInput`）。
- **四面体 × 输入**：BVH 剪枝 + 顶点包含 / 棱-三角形 / 三角形-棱 四类测试（`tetIntersectsInput`）。
- **偏移面**：完全隐式，只有 `d(x) − δ` 的求根，从不显式构造网格。

## 4. 验证结果

用 `test/validate_obj.py` 检查：边关联数、朝向一致性、有符号体积、输出顶点到输入距离、以及输入顶点相对输出网格的绕数（≈1 表示被严格包住）。

| 模型 | 输入三角面 | α / δ | 输出 | 时间 | 边关联 | 顶点距离 | 绕数 |
|---|---|---|---|---|---|---|---|
| 球（icosphere） | 1,280 | 0.04 / 0.004 | 6.2k V / 12.5k F | 6.2 s | 全 2 | 全 = 0.004 | 1.000 |
| 环面 | 2,880 | 0.05 / 0.004 | 4.6k V / 9.2k F | 19 s | 全 2 | 全 = 0.004 | 1.000 |
| 盒子（尖边/共面） | 12 | 0.1 / 0.01 | 1.8k F | 0.06 s | 全 2 | 全 = 0.01 | 1.000 |
| 开放平面（双侧包裹） | 800 | 0.1 / 0.02 | 2.4k F | 0.6 s | 全 2 | 全 = 0.02 | 1.000 |
| 缺陷球汤（挖洞+翻面+扰动） | 1,220 | 0.04 / 0.004 | 21.8k V / 43.9k F | 50 s | 2/4 | 全 = 0.004 | 1.000 |

说明：

- 球/环/盒的“边关联”全为 2，即闭合 2-流形；缺陷球汤出现 **21 条 4-关联边**，是论文明确指出的 *immersed manifold* 的“pinching”，不是破洞（没有奇数关联边）。
- 开放平面输出体积 ≈ `面积 × 2δ`，即论文的 two-sided wrap。
- 盒子的输出体积略大于输入体积，且所有顶点距离严格等于 δ。

## 5. 限制

1. **非精确谓词**：`long double` + 相对过滤，不是 Shewchuk 自适应精确谓词。极端共球/共面输入可能被插入回退的微小抖动（≤ `1e-8·diag`）扰动；输出不保证与理论完全一致。
2. **性能**：单线程、无空间排序、BVH 是朴素中位数分割。比 CGAL 的 `Alpha_wrap_3` 慢一个量级左右；10 万级输出面片可能需要分钟级。
3. **内存**：四面体只追加不回收（死亡单元保留），大模型内存会持续增长。
4. **参数全局**：α、δ 是常数，没有自适应 α 场。
5. **pinching**：与论文一致，输出是 immersed 组合流形，允许边处自接触；若要严格嵌入需在雕刻时禁止该类配置（本实现未做）。
6. **输入规模**：BVH 与 Delaunay 都是单线程；百万级输入建议先简化或分块。

## 6. 调参建议

- 先 `α ≈ diag/50`、`δ ≈ diag/500` 试跑；
- 要更贴合输入：减小 δ（注意 R2 变多、质量下降、时间上升）；
- 要更粗/更快：同时增大 α 与 δ；
- 要保留薄壁、缝隙、小孔：α 调到小于缝隙宽度；
- 输出出现“包住内部空腔”的双层结构是正常现象（α 小于开口尺度时，front 会钻进孔洞）。

## 7. 目录

```
src/alpha_wrap.cpp      单文件实现
CMakeLists.txt
build.bat
test/gen_models.py      生成球/环/盒/平面/缺陷汤（需要 numpy）
test/validate_obj.py    输出网格校验（水密、朝向、体积、δ、绕数）
test/inspect_obj.py     统计 OBJ 顶点/面数/包围盒
test/render_obj.py      软件光栅化预览图（需要 Pillow）
test/*.obj              测试输入（生成的 *_wrap*.obj 由 .gitignore 排除）
```

## 8. 参考

C. Portaneri, M. Rouxel-Labbé, M. Hemmer, D. Cohen-Steiner, P. Alliez.
*Alpha Wrapping with an Offset.* ACM Transactions on Graphics 41(4), Article 127, 2022.
