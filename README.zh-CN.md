# Reinhardt 最大周长小多边形：n=16、32、64 计算机辅助证明候选

本仓库包含 Reinhardt 最大周长小多边形问题在

\[
n=16,\quad n=32,\quad n=64
\]

三个二次幂情形中的证明候选、可执行验证程序和交叉复算材料。

> **状态（2026 年 8 月 5 日）：**这些仍是证明候选，不是已经经过同行评审确认的定理。证明开发过程重度依赖 **OpenAI GPT-5.6-Sol**；它在 Jizhou Guo 的反复提示和反馈下生成了数学论证与验证代码的相当大部分。三套原始证书均已复现，后续 AI 辅助交叉检查目前没有发现具体错误，但仍需要独立的人类专家审查。

## 当前结果

| 情形 | 候选最优周长 | 全枚举覆盖 | 幸存编码 |
|---|---:|---:|---:|
| \(n=16\) | `3.136547716486607386085967...` | 全部 \(2^{15}=32768\) 个规范化编码 | 16 个，一个二面体轨道 |
| \(n=32\) | `3.1403311569546193658254013805774586723...` | 全部 \(2^{31}\) 个规范化编码 | 96 个，三个轨道 |
| \(n=64\) | `3.1412772509327728680619914155024682980...` | 通过严格三元 meet-in-the-middle 覆盖全部 \(2^{64}\) 个半编码 | 896 个，六个轨道 |

三份证明候选均声称最优合同类唯一。

## 运行验证程序

环境要求：Python 3.10 以上；完整 \(n=64\) 扫描需要 C++17 编译器和约 2 GB 可用内存。

快速验证：

```bash
make verify-fast
```

完整验证，包括重新生成 \(n=64\) 的 896 个幸存编码：

```bash
make verify-all
```

核心证书使用严格的有理数/整数区间判定。`audits/` 中的程序属于不同实现和高精度交叉复算。

## 文件结构

```text
cases/      证明候选、verifier 和记录输出
audits/     独立扫描和高精度复算
docs/       状态、审查重点、AI 披露和参考文献
scripts/    一键复现脚本
```

## 当前最需要审查的部分

计算证书已经得到较强复现；真正需要领域专家逐行检查的是：差体重构与饱和性、近正则局部化、离散谱估计，以及固定编码下 KKT 点的唯一性。详见 [docs/REVIEW_REQUEST.md](docs/REVIEW_REQUEST.md)。

## AI 使用与人类贡献

这项工作并非只是用 AI 润色。证明开发重度依赖 GPT-5.6-Sol。Jizhou Guo 主要负责选择问题、提供 prompt 与资料、反复引导、整理生成材料并发起后续审计；不声称独立推导了每一个证明步骤或亲手编写了每一行代码。后续复核还使用了 Claude 和 GPT-5.6 Thinking，但这些仍不等于独立人类同行评审。完整说明见 [docs/AUTHORSHIP_AND_AI_DISCLOSURE.md](docs/AUTHORSHIP_AND_AI_DISCLOSURE.md)。

## 维护者

**Jizhou Guo**

- [ORCID](https://orcid.org/0009-0001-0699-9164)
- [Google Scholar](https://scholar.google.com/citations?user=fcBDdsYAAAAJ)
- [DBLP](https://dblp.org/pid/378/4049.html)
- [个人主页](https://aster2024.github.io/)
- [X / Twitter](https://twitter.com/TheOsmanthus)

## 许可证

全部内容采用 MIT License，详见 [LICENSE](LICENSE)。
