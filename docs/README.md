# 分析記錄索引

每一次有結論的分析都存成一份 md，按日期排序。
**如果後續的開發或實驗改變了某個結論，回去更新對應的檔案，並在檔案開頭的「狀態」欄註明。**

| 日期 | 檔案 | 主題 | 狀態 |
|---|---|---|---|
| 2026-08-15 | [01-camera-characterisation.md](01-camera-characterisation.md) | 從影像反推相機參數：tone curve、自動白平衡、自動曝光 | 有效 |
| 2026-08-15 | [02-all-fields-confound.md](02-all-fields-confound.md) | 12 個 field 全跑 + 混淆因子分析 | **有效，且推翻了先前的 T/C = 1.66** |

## 目前的總結論（一句話）

> 這批 HMC3 / MyD88 影像**無法區分「treatment 效應」與「細胞密度經由自動曝光造成的假象」**。
> 在補到必要的對照資料之前，不應該報告任何效應量。

詳見 [02-all-fields-confound.md](02-all-fields-confound.md)。

## 待補資料（依重要性排序）

1. **control / treatment 的分組** — 目前完全未知，先前沿用的假設已被時間戳推翻
2. **secondary-only（只打二抗）對照影像** — 定出加性 offset
3. **相機當時的實際設定** — 使用者將到實驗室確認後回報
4. **相機 tone curve 校正 bracket** — M 模式、固定場景、已知曝光階梯
