# SAP v0.1（HL 修訂版）檢視紀錄

檢視日期：2026-10-01
檢視對象：`APOE × CKD Progression — Statistical Analysis Plan ..._HL.docx`（與線上文件 rev 136 內容一致）

HL 版將分析範圍收斂為 JAMA 2005 的主體複製：刪除 summary score、genetic PC、親緣處理、biobank 延遲進入、multiple imputation、Model 4、次族群分析，以及 eGFR trajectory／Dementia／Kremezin（原第 8 節）。程式碼依此範圍撰寫。以下為仍需決定或補寫之處，依影響程度排序；「程式預設」欄為程式碼目前的處理方式，可於 `config/config.yml` 修改。

## A. 影響主要結果、建議補入 SAP

| # | 位置 | 問題 | 建議 | 程式預設 |
| --- | --- | --- | --- | --- |
| A1 | §1、§6.3 | 刪除「主要檢定（1 項）」後，主要檢定未明確指定。§6 寫「2-df 檢定」，§6.6 檢定力卻以「每 ε4 對偶 HR」計算，兩者不一致。 | 指定唯一主要檢定。建議：Model 3 之 ε2 + ε4 2-df LRT 為主要檢定，ε2、ε4 每對偶 HR 為主要效應量；檢定力改以此為準。 | 兩者皆輸出；Table 2 以 2-df LRT P 值呈現 |
| A2 | §2 | 移除親緣處理。iHi 親緣比例約 30–40%，未處理會低估標準誤，且 APOE 在家族內共享。 | 至少擇一：(a) 每家族保留一人；(b) 以家族為 cluster 之 robust SE。 | `kinship.method: none`；可改為 `cluster` 或 `unrelated`，需提供 `family_id` |
| A3 | §2 | 移除 biobank 延遲進入。若 index 早於 iHi 收案（同意書）日，參與者必須存活至收案日才能被納入，形成 immortal time；ε4 若與死亡相關，偏誤方向無法預期。 | 至少列為敏感度分析：以 max(index, 同意書日) 為進入時點（left truncation）。 | `delayed_entry: false`；設 `true` 並提供 `consent_date` 即啟用 |
| A4 | §5、§6 | 刪除 imputation 與 complete-case 敏感度分析後，SAP 未說明遺漏值處理。BMI、血脂、血壓在 EMR 中遺漏比例可能偏高。 | 明定 complete-case 為主要分析，並報告各共變項遺漏比例；若 Model 3 遺漏 >10%，再考慮 MI。 | Complete-case；Model 1–3 共用 Model 3 完整樣本（`common_sample: true`），另輸出遺漏比例表 |
| A5 | §4 | 「死亡資料」段落已清空，但主要 outcome 第 2 項與設限仍使用死亡。EMR 僅有院內死亡時，院外死亡被視為失聯設限。 | 確認可否串接衛福部死因檔；若無，明寫「死亡 = 院內死亡」並列為限制。 | `person.death_date`；CKD 死亡 = 死因碼（若有）或死亡前 30 天內有 CKD 診斷 |
| A6 | §4 | Baseline 已有 CKD 診斷碼者，追蹤期間任何住院只要帶 N18 次診斷即成為事件，主要 outcome 會混入 prevalent CKD。ARIC 收案時 CKD 住院史極少，此問題在 EMR 中較嚴重。 | 擇一：(a) 排除 index 前已有 CKD outcome 碼者；(b) 保留但做排除之敏感度分析。 | 保留；另跑敏感度分析「排除 index 前 CKD 診斷」 |
| A7 | §4 | 肌酸酐上升之「確認」規則不完整：(i) 確認值是「≥90 天後下一筆」還是「≥90 天後任一筆」；(ii) 候選事件後無後續檢驗（死亡、失聯）時是否算事件。 | 建議：下一筆；無後續值者不算事件，另以「不需確認」敏感度分析涵蓋。 | `confirm_rule: next`；無後續值 = 非事件 |

## B. 定義需補細節

| # | 位置 | 問題 | 建議 |
| --- | --- | --- | --- |
| B1 | §3 表格 | 「稀有型別(ε1)估計不穩」：兩個 SNP 組成的六種基因型不含 ε1。ε1 會在雙異合子（rs429358 TC、rs7412 CT）中與 ε2/ε4 混淆，未 phasing 時通常歸為 ε2/ε4；其餘組合（如 rs429358 CC + rs7412 CT）才確定含 ε1。 | 改寫為：「稀有基因型（ε2/ε2、ε4/ε4）估計不穩；含 ε1 之組合排除，雙異合子依慣例歸為 ε2/ε4」。程式依此處理並輸出 ε1 人數。 |
| B2 | §4 | CKD 住院事件日未定義。JAMA 2005 使用出院日。 | 明定出院日（若無則入院日）。程式預設 `hosp_event_date: discharge`。 |
| B3 | §4 次要 outcome | 次要 outcome 表保留，但 §6 已無分析方式（原 7.6 已刪）。 | 補一句：「次要 outcome 以 Model 3 對偶基因加成模型分析」。程式依此執行。 |
| B4 | §4 次要 outcome | ESKD「長期透析 ≥90 天」之操作定義未明。 | 建議：首次透析處置且 ≥90 天後仍有透析紀錄，或腎移植；重大傷病（ESRD）如可串接則優先。 |
| B5 | §5 共變項 | BMI 未定窗口；糖尿病之檢驗與用藥窗口未明。 | 程式採：BMI 同血壓 ±180 天；HbA1c／血糖同檢驗窗口（−365 至 +30 天）；降血糖藥（ATC A10）同用藥窗口（前 90 天）。 |
| B6 | §6.3 | 「違反比例風險時以時間分段估計」未指定切點。 | 程式預設切點 3、6 年（可於 config 修改）。 |
| B7 | §6.5 | 「年齡為時間尺度」之 left truncation 需以 index 年齡為進入時點，模型移除年齡共變項。 | 程式依此執行。 |
| B8 | §6.6 | 檢定力計算核對：p = 0.09、HR = 0.85、α = 0.05、power 80% → 1,814 件事件，與 SAP「約 1,800」一致。此計算未考慮與 ε2 共同入模，影響可忽略。 | 無需修改；盤點事件數後以 `11_power.R` 反推可偵測 HR。 |

## C. 範圍與一致性

- 原第 8 節（eGFR trajectory、Dementia／CVD、Kremezin effect modification）已刪除，但 9/14 會議與 9/15 任務指派之主要問題即為 eGFR trajectory 與 Kremezin。建議另立 SAP（或附錄），避免日後被視為事後分析。本 repo 未納入這三項。
- 刪除 genetic PC 可接受（iHi 以漢人為主），但若投稿遺傳流行病學期刊，審稿人可能要求，建議保留 PC1–10 為敏感度分析。程式中為選項 `adjust_pcs`，預設關閉。
- §7 待決事項中「CVD 是否納入」在 HL 版已無對應分析，可刪除或移至新 SAP。
