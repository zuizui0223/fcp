# FCP / RGFCA — 12種・実写真1,034枚の観察分布と花色ROI技術監査（2026-10-16）

**種内花色多型の正式な生態学的結果ではない。** この12種は 2026-09-27以前のcreated-atフィルタを持つlive iNaturalist APIから取得した技術コホート。写真IDは厳密な2026-09-27月次bulk Open Datasetアーカイブと照合されていない。原位置精度・観察日もこの技術データに含まれていない。したがって、以下は写真を通じて記録された観察地理であって原産／自然分布域ではない。

## 本当に実行・検証した範囲

- Original ROI-v4/EfficientSAM実画素測定: [wave1 #38832829706](https://github.com/zuizui0223/fcp/actions/runs/38832829706)（552枚・分類可能202枚）と[wave2 #38840579008](https://github.com/zuizui0223/fcp/actions/runs/38840579008)（482枚・分類可能203枚）。合計 **1,034枚／12種・分類可能405枚**。
- 元photo ID・observation ID・observer ID・sourceの座標・元ROI terminalのSHAを照合する独立再計算：[Actions #38862008560](https://github.com/zuizui0223/fcp/actions/runs/38862008560) **SUCCESS**、artifact `12030290004`。今回画像を再ダウンロードしておらず、古い花色ラベルも利用していない。
- 派生12種CSV SHA256: `705b0f13f32c0fbfcbb1c811efbc2dc9f1313883ffe96d1bfef26b0147fccbeb`。再現コード: `scripts/analysis/audit_rgfca_12_original_live_photo_geographic_span.py`。

## 元画像観察の実際の地理分散と花色品質

同じ観察IDに複数写真があっても1観察として数えた。独立観察3件以上・観察者2人以上の5°セルの重心を球面大円距離で比較した。距離は**位置精度未検証の技術値**である。1,000km以上の広域候補は4種、95%重心距離200km以内の狭域候補は2種、残り6種は中間・未解決。元写真の枚数上限と地理観察量があるため、真の分布幅への外挿は不可。

| 種 | 支持セル間最大距離 km | ROI分類可能/実写真 | 同一ROI内に白・黄各10%以上 | 暫定地理 |
|---|---:|---:|---:|---|
| *Pseudodictamnus undulatus* | 155 | 25/96 | 6 | 狭域候補 |
| *Acourtia nana* | 1,351 | 31/96 | 2 | 広域候補 |
| *Eriogonum villiflorum* | 162 | 22/70 | 0 | 中間・不明 |
| *Thelesperma filifolium* | 630 | 67/96 | 7 | 中間・不明 |
| *Fimbristylis alboviridis* | 0 | 6/30 | 1 | 中間・不明 |
| *Combretum grandidieri* | 191 | 31/89 | 0 | 狭域候補 |
| *Chloris formosana* | 337 | 19/96 | 2 | 中間・不明 |
| *Lactuca perennis* | 1,126 | 40/96 | 2 | 広域候補 |
| *Matricaria chamomilla* | 19,412 | 64/96 | 58 | 広域候補 |
| *Azorella hydrocotyloides* | 228 | 43/94 | 17 | 中間・不明 |
| *Nassauvia axillaris* | 2,196 | 43/96 | 1 | 広域候補 |
| *Grevillea decurrens* | 467 | 14/79 | 1 | 中間・不明 |

## 解釈と使用上の禁止事項

- *Matricaria chamomilla* は64枚の分類可能な花ROI中58枚で白・黄が同一ROIに共存していた。頭状花序の花器官による固定的な配色などを除外せず、「色ラベルが反復して現れたから遺伝的多型」と解釈してはいけない。元写真クロップの人手確認は未実施。
- *Lactuca perennis* は40枚の分類可能写真に白30・青紫6・黄橙3・赤桃1という複数ラベルが観察されているが、撮影者・季節・花器官・個体の独立性を満たす局所間差と断定できない。
- 観察地域が遠い写真の存在と『複数地域で花色型が入れ替わった』ことは別々の検定である。正式な広域／狭域源コホート1,000種の分類、P_poly、W/H2、気候／土壌の生態効果の値をここから作らない。

**次の主線:** [真正2026-09-27月次原本の1,000種個別photo ID凍結 #38847189688](https://github.com/zuizui0223/fcp/actions/runs/38847189688) が成功した場合にのみ、別途コード化済みの位置精度<=50km・観察者独立の地理分散→狭域対照ペア→元ROI画素測定の正式実行へ移す。技術12種の画像IDと花色terminalを本番には再利用しない。
