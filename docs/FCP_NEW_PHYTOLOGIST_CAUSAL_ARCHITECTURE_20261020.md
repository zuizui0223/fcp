# FCP / New Phytologist：世界的花色変異の「生起・配置・維持」を区別する仮説体系（2026-10-20）

**状態：研究上の再構成案。** 既存の `docs/POLYMORPHISM_MANUSCRIPT_NEW_PHYTOLOGIST.md` は変更しない。原稿改訂方針と、既存447,500枚実測データから支持されること／新規解析が必要なことを明示する。

## 中心問題と提案する理論

> **なぜ被子植物の異なる種で似た花色の分化方向が繰り返し現れるのに、その地理的な配置と維持機構は異なるのか？**
>
> Why do recurrent directions of intraspecific flower-colour variation yield different geographic mosaics and persistence across plant species?

「白が暖地で有利」でも「花色が空間構造を持つ」でもなく、**表現型の生起（developmental accessibility）→ 地域間の配置（historical/demographic + environmental filtering）→ 集団内の維持（fitness, pleiotropy, mating）**の3段階を比較する仮説。第2段階でhistoryがenvironmentより主要かは、現時点では**未確定**。「IBD-like > IBE-like」は分布史そのものの測定ではない。

## 1. 結論として許される確定／探索的結果（混ぜない）

| 問い | 保存済み結果 | 状態／上限 |
|---|---|---|
| 実測画像 | 4,475種×各100枚＝447,500枚、分類可能227,107枚 | 同種写真の色型ラベル。野生の遺伝型を未確認 |
| 色型空間の共通方向 | 凍結された白―有色方向 `W=0.517` 対構成保存帰無中央値 `0.457`、`p=.001`（前向き種・写真非重複） | 同一撮影系内の事前固定H2結果。ただし白ラベルは白飛びと関連し、厳格な技術妥当性ゲートは未確定 |
| 地域的色構造 | 50 km: 1,761種の平均 `D_global−D_local=+0.017615`; 正の効果は984/1761（55.88%）、対応する固定組成帰無平均44.58% | 元々の色頻度・地点を固定した写真ラベル置換より地域内の色差が小さい。個別984種の検定ではない |
| 白―有色と有色―有色 | 地理構造の83.81%は白―有色由来、だが全域異色ペアの79.57%が元々白―有色。規格化（両成分baseline>1%の636種）で白+0.05696、有色間+0.06945、差95%CIに0 | **白―有色が特別に強く地理分化するとは支持されない**。種の97.8%に白ラベルがあることは色分類・画像効果から独立ではない |
| 気温 | BIO5と白の種内関連: 1,574種・86,114写真、1℃当たり種等重み平均log-odds `+0.03233`、属bootstrap95% `[+0.01786,+0.04728]` | 小さい正の**位置平年値との関連**、post hoc。古い事前固定コホート間転移条件は非達成。適応的白化の証明ではない |
| IBD-like と IBE-like | 原則の3コホート：geography | BIO5 差を調整後の種平均partial rho `0.02760/0.02414/0.02225`、BIO5 | geography のpartial rho `0.00925/0.00826/0.01019` | **IBD-like は3群で IBE-like より大きい**。ただし表現型の距離相関であり、真の遺伝的IBD/IBEでも、歴史の寄与率でもない。BIO5残差は花−背景・厳格観察者検証では弱い |
| 静的分布域・ニッチの要約 | 別の文献由来C/S比較（v6）で、range size・GBIF fragmentation・単純静的環境turnover・climate breadthは頑健な説明因子にならない | 別データ、44.75万枚の結果と混同不可。**歴史優位の実証ではない** |
| 同一局所集団でのC共存 | 同一月・別観察者・5 km条件の帰無テストで支持なし（広い条件で194種 p=0.300、厳条件14種 p=0.275） | 50 kmの均一性 ≠ 遺伝的FCP共存／固定 |
| 分子／維持 | *Silene littorea* のPAL 8–21%、WAL <1%（原著の21集団）；異色RNA-seq候補・器官別色素など | **公開個別研究と文献表の証拠**、FCP画像個体に直接結合した遺伝・分子試験ではない |

## 2. 論文の五つの競合／階層仮説

### H1 世界の形質「生起の偏り」：Constrained phenotype accessibility

**問い**：種ごとに分子背景が違っても、白―有色という変化が他の色相転換より共通しやすいのか。

- 予測 H1a：既存4色ラベルの頻度や白飛びを保存した適切な構築帰無より、連続パレットの白―有色方向に過剰な種内変化が並ぶ。
- 反証：白い花に見えるのは撮影露出・花以外の白い部分・粗い分類のためだけであり、花器官の独立再測色では過剰方向が消える。
- 既存の直接支持：元のprospective W 検定。**新規課題は**独立画像系で器官限定・色彩校正した再検証。

### H2 歴史・移動のフィルタ：Historical sorting / dispersal limitation

**問い**：地理的な花色分化の多くが、同一種内での分布域変遷と接続性の違いによって生じるのか。

- 予測 H2a：環境距離が同程度でも、独立した地史的分断、障壁、推定中立遺伝構造の隔たりが大きい集団ペアほど色差が大きい。
- 予測 H2b：拡大前線と氷期レフュージア由来集団など、過去の分布変遷が異なるセクター間では、現在の同じ気候下でも色型組成が異なる。
- 反証：独立の中立集団構造／歴史指標を制御した後に距離効果が消え、環境不一致が主説明力を持つ。また、遠いが歴史的に接続した地域では色型組成が一致する。
- 既存結果：表現型IBD-like はIBE-likeより大きいが、これは**歴史を直接検証していない**。既存のGBIF fragmentation proxy は別の小標本で不支持。
- 必須データ：標本位置・個体群の中立SNP構造、系統地理、歴史・LGM分布変遷は独立指標として候補化。LGMのSDM後方投影は不確実なproxyであり、コア歴史証拠と同一視しない。

### H3 現在環境のフィルタ：Environmental sorting conditional on geography/history

**問い**：地理的に隔離された同種集団の色差は、気候・土壌・送粉／食害環境によって独立に組織化されるのか。

- 予測 H3a：地理距離・調査機会・季節・観察者を合わせた同種集団ペアでも、環境差が大きいほど花色差が大きい（phenotypic IBE-like）。逆に同じ環境にある遠方集団の色は相似になり得る。
- 予測 H3b：BIO5の白方向正の平均関連を保持するが、乾燥・土壌・訪花ギルド等との交互作用で種ごとの符号・強さが変わる。**すべての種が暖地で白化するとは予測しない**。
- 反証：同等地理距離と観察条件で環境差を変えても色差が残らない、または中立構造で説明され尽くす。
- 既存結果：post hoc BIO5 `+0.03233/°C`、弱い部分順位 IBE-like。元の厳格背景／観察者判定は未支持なので強い因果主張不可。

### H4 空間と時間の対比：Coexistence versus spatial segregation

**問い**：同じ遺伝的色型が局所集団で維持される場合（C）と地域間に分かれる場合（S）とで、何が違うのか。

- 予測 H4a：**持続的な集団間の環境差＋遺伝子流動の制限**でSが増え、**局所的に時間変動する選択・複数送粉者・自殖保障／種子バンク**ではCが維持され得る。双方が作用するとC+Sのモザイクになる。
- 反証：局所複数年・個体／遺伝型確認のあるC/S分類でも、これらの環境・歴史説明変数と無関係。
- 重要：現在の50 km統計からC/Sを直接分類しない。厳格な同一局所集団での撮影・個体追跡・遺伝型、ERA5-Landなど年別気候と送粉者時系列が必要。BIO4は年間季節性であり年度間変動と異なる。

### H5 多型の維持に伴う器官特異的コスト：Tissue-specific pigment control filters persistence

**問い**：白花が「生じる」ことと「野生集団で維持される」ことを隔てるのは、花弁だけの色素制御か、それとも葉・茎に及ぶ負の多面発現か。

- 予測 H5a：花弁のアントシアニンを低下させ、葉・茎の防御色素を保持する**PAL**は、葉まで失う**WAL**より高頻度で長期維持され得る。
- 予測 H5b：PAL/WALの頻度差は送粉者・乾燥・食害・自殖率・個体群史で修飾される。白花そのものを単一の色素代謝コスト・利益へ還元しない。
- 反証：独立遺伝背景・同じ環境で比較したPAL/WALの生存・繁殖差や長期集団頻度が維持予測と対応しない。
- **直接の公開機構例**：
  - *Parrya nudicaulis*: 白花弁で`CHS`発現約24倍低下、葉・萼は同等（petal-limited）；白の頻度と暖季気温に正の関係。**因果DNA変異は未同定**。doi:10.1371/journal.pone.0018230。
  - *Silene littorea*: 白花蕾`F3h1`発現は有色より >42倍低く、`Myb1a`が候補；拡張SNP調査で色型を一貫して決める変異は未発見。doi:10.3389/fpls.2016.00204。
  - 同種 *S. littorea*：野生21集団でPAL 8–21%（存在する2集団）、WAL <1%。器官別HPLCでPAL/WALを識別。doi:10.1186/s12870-019-2082-6。
  - *Orychophragmus violaceus*：F2紫929／白295が単一劣性座位と整合、1.37 Mb範囲で`OvANS`候補。ただし候補遺伝子の因果機能証明ではない。doi:10.3390/biology14121669。
  - **反例が重要**：*Boechera stricta*では乾燥・低標高で着色が誘導され、十分な水条件では白花系統に繁殖上の利益が認められた。doi:10.1111/nph.14998。したがって全球のBIO5正平均を全種共通選択方向へと読み替えられない。
- **研究の直接の未解決点**：これら論文の花弁／葉の分子データと、44.75万写真の個体／遺伝型は直接連結していない。集計したPAL13／WAL13文献表も無作為標本ではない。写真の白ラベルが必ずアントシアニン欠損であるわけではない。

## 3. IBD / IBEを「歴史 vs 環境」の検定にする条件

既存の研究は、連続花色距離と（a）地理距離、（b）BIO5距離の偏順位相関を比較した。

| 元の三コホート | phenotypic IBD-like `rho_geo|BIO5` | phenotypic IBE-like `rho_BIO5|geo` |
|---|---:|---:|
| Discovery | +0.02760 | +0.00925 |
| Validation | +0.02414 | +0.00826 |
| Prospective | +0.02225 | +0.01019 |

両成分は正だが小さく、地理の方が大きい。**「分布史が環境の約3倍を説明した」という解釈は不可。** 地理は時空間に沿った未測定環境、撮影努力、分散障壁も吸収する。元のIBE-likeは花／背景・同一観察者検証で不支持だった。既存の結果を保持したまま、これを次の対照へ発展させる：

1. **写真由来の表現型IBD/IBE**を第4・第5コホートまで拡張し、花器官・露出・栽培／野生・花期・季節を統制する。結果の主単位は種等重み、属・科層と補正付きCI。
2. **同じ地理距離の集団ペア**で環境差が小さい／大きいペアを比較。加えて、地理的に遠いが類似環境のペアを対照に取る（背景環境availabilityが異なるときは補正）。
3. **独立の中立遺伝距離 `F_ST`／SNPクラスタ／現生接続性**を入れて `colour ~ neutral genetic distance + environmental distance + geography + sample structure`。系統内で集団背景と花色座位を明示的に分ける。
4. **分布域変遷・レフュージア／拡大方向**は独立の地史データ・過去分布モデルで感度付きに評価。現存range areaやGBIF fragmentationが非有意だったことを「歴史は重要でない」とはしない。
5. **条件付け後の残差・機構の方向**を判断。IBD/IBEの数値だけでどちらが因果主導か結論しない。

## 4. Global patternを最初に見せるには何が不足しているか

「全球4,475種」というデータ規模と「全球の生物地理区によって多型発生率が違う」という結論は**別物**。現状は4,475種が各100枚の機会標本で、栽培／外来記録・白飛び・地域／観察者偏りを残す。したがって完成したのは種内**写真色分布の比較**であって、亜熱帯・山地・島嶼の多型が多いという全球生物地理パターンではない。

最初の実証図を作る優先操作：
- 各種の撮影地点をbiogeographic realm、climate zone、島嶼／大陸、標高帯で層別し、**種−地域当たりの有効写真数・観察者数を維持**して現地色ラベル多様性・標準化Dを推定。外来／栽培候補と確定野生は分離。
- 「多型率」ではなく**確認された同種地域内花色の共存／地理的置換／未確定**という証拠水準つき世界地図を作る。5 km同一花期独立観察者条件は現状非支持で、Cを無理に自動判定しない。
- **適応推論ではなく空間的形質差の強さ**と、地域別の方向、温度・距離との独立な対応を可視化。セル別種等重み・地理機会null、同数レアファクション、残存カラー分類誤差を使う。
- 世界共通の境界hotspot集中も検討可能だが、既存RGFCAの集中主検定はp≈0.07で非支持だった。肯定的世界地図を後付けで選ばない。
- 系統推論：種内の比較で種間色の基本差は統制されるが、各種の傾きの一般化では共有祖先を再考慮する（属・科クラスター／独立系統樹）。

## 5. 論文のFigと原稿の順序（測れていること／未測定を明記）

**Proposed title, hypothesis-led:** *Recurrent floral colour variants, divergent geographic histories: an integrative framework for intraspecific flower-colour variation*。歴史を直接実測するまで、このタイトルの後半は**仮説**であり既存結果だけで投稿する安全なタイトルにはしない。

1. **Figure 1 | Global photographic phenotype variation atlas** — 地域・分類群での種内色多様性と白―有色軸の頻度（地域層別・努力補正・native filterが追加実行必要）。現状は全種の分母とWがあるのみ。
2. **Figure 2 | Geographic allocation rather than mere presence** — 50 kmの種別Δとmatched null／地理的clines・mosaics。5 kmのC判定は機能的に未確認と明記。**既存実測の主結果**。
3. **Figure 3 | Historical-distance and environmental sorting as competing explanations** — IBD-like vs IBE-like（3コホート実測あり）、5コホートへ拡張し独立地史/中立遺伝距離の層を追加できればmechanistic claimに進める。**現状の歴史優位は未証明**。
4. **Figure 4 | Conditional climate association and ecological sign heterogeneity** — 1,574種の平均BIO5正傾向と各種分布、環境×ギルド／環境幅を測る（平均は実測、関係の異質性・交互作用は未確定）。
5. **Figure 5 | What makes a white morph persistent?** — PAL/WAL organ-specific pigment mechanism, case evidence anchored in *Parrya*, *Silene*, *Orychophragmus*, *Boechera*. 現状は**別研究の分子・生化学の統合**であり44.75万枚由来の分子実測ではない。因果分子結合の主論文図として昇格するには同じ自然集団・個体のgenotype、色素量とfitnessが必要。

**現実的な出版線引き：**
- **今の44.75万枚のみで言えるNPの結果：** 再現する連続白―有色軸 + 多数種の地理的配置 + 小さなBIO5正平均 + phenotypic IBD-like優位（最初の3群のみ） + 異なる維持機構の候補。
- **より強いNPの因果主張に足りないもの：** 真正なregional-scale world map、独立の分布史／中立遺伝構造、5コホートのIBD/IBE厳格空間条件比較、器官別色素と適応度。
- 先行文献から機構を提案することは可能だが、新しい分子データの実験・遺伝子型なしに「分布史とPALのコストが多型維持を決める」とは言わない。

## 6. 実行順（主線）

1. **47万ではなく44.75万枚を使った全種の大陸・biome・native/wild品質別の地理atlas**：種−地域機会分母と色ラベル品質を固定し、独立画像検証可能な色差だけ地図化する。
2. **第4・第5コホートへの連続色IBD/IBE拡張と同距離・異環境ペア比較**：3群での地理優位が広い種群でも再現するか、環境差が独立な説明力を持つかを確認。
3. **独立の歴史・接続性 proxy の取得**：この段階で初めてhistory > current environmentの仮説を検定する。単純なGBIF断片度は既に別標本で不支持なので再利用だけで終えない。
4. **PAL/WALを軸にした複数種の分子・器官別証拠の統合**：同種画像個体との一致、器官、色素経路、因果アレルか発現相関か、fitness測定の有無を階層化。遺伝的分化の機構が明瞭な少数系を現地検証へ移す。

## 根拠・参照元

- 原論文の凍結表現型方向・技術問題：`docs/POLYMORPHISM_MANUSCRIPT_NEW_PHYTOLOGIST.md` と `docs/POLYMORPHISM_SECONDARY_MECHANISM_EVIDENCE_LEDGER_20260928.md`。
- 44.75万枚・1,761種：`docs/FCP_447500_SPECIES_DISTRIBUTION_20261018.md`、`results/fcp_real_447500_white_axis_decomposition_20261018/component_summary.csv`。
- 1,574種BIO5：`docs/FCP_COMMON_UNITS_TRANSPORT_MOLECULAR_20261019.md`。
- 元の3群IBD/IBE：`results/polymorphism_phenotypic_IBD_IBE_posthoc_20261007/result.json`。
- 787種の温度傾き強度×空間構造の探索診断：`docs/FCP_TEMPERATURE_MAGNITUDE_AND_GEOGRAPHIC_SORTING_20261019.md`。絶対値と地理構造の偏相関は約+0.0914だが同じ元写真由来であり共同ラベル置換なしの独立機構証拠ではない。
- C/S別サンプルの静的proxy非支持：`docs/JBI_FCP_MECHANISM_HYPOTHESES.md`、`docs/JBI_SPACE_TIME_HYPOTHESES.md`。
- *Boechera*: Vaidya et al. (2018), New Phytologist, https://doi.org/10.1111/nph.14998
- *Parrya*: Dick et al. (2011), PLOS ONE, https://doi.org/10.1371/journal.pone.0018230
- *Silene* F3H/MYB: Casimiro-Soriguer et al. (2016), https://doi.org/10.3389/fpls.2016.00204
- *Silene* PAL/WAL: Del Valle et al. (2019), https://doi.org/10.1186/s12870-019-2082-6
- *Orychophragmus*: Li et al. (2025), https://doi.org/10.3390/biology14121669

**Original frozen prospective confirmations remain separate. No new environmental data, genotypes, image classifiers, historical distribution reconstructions or causal experiments were generated to produce this architectural note.**
