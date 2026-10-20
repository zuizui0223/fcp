# FCP H1–H5：最小追加検証の判別的プロトコル

**Freeze date:** 2026-10-20 (JST)  
**Research question:** *Why do flowering plants repeatedly diversify along similar colour axes, yet differ in how and where that variation is maintained?*  
**Scope:** intraspecific flower-colour variation; photographs ≠ genetically inherited sympatric colour polymorphism.  
**Original main / 158-species prospective W / 447,500 photos / 1,761-species 50-km null are immutable.**

## 結論と追加の最小単位

既存の447,500枚から世界の写真花色・地域分布を見せる。新規の追加作業は以下の**二種類だけ**に集約する。

**A. H1に対する独立測色**：既存W選定規則を凍結し、写真以外の観察者／センサー／校正方法による独立花器官色計測で、Wが単なる撮影・ラベル構築由来か試す。対象候補20種（少なくとも10科が目標）、各種2色型×8個体以上を目標とするが、**新データを開く前に対象種と欠測方針を固定**。既存の花色結果の強弱で都合の良い種を選び直さない。集めるのは個体単位の主花弁反射率／色校正写真・生育地点／花齢／光条件・標本由来。選択後の20種でWが支持されても、全球で白色多型が何%かは推定できない。

**B. H2–H5を共通の野生集団で試す連結パイロット**：生物学的な白／有色型が現地で確認され、採取許可と繁殖型が確認できる**独立2種**を対象とする。目安は**各種6個体群×各個体群20個体＝240標識個体**、少なくとも2開花季節、集団間で現在の生育環境の似た組合せ／違う組合せを採る。個体群内に二色型が確認できない場合でも、その集団を恣意的に排除しない。空間選択の比較と遺伝構造の識別は**集団／独立な家系**を単位とし、写真ペアや同一個体の花を独立反復として数えない。年2回採取は2年間の再観測に過ぎず、何世代も続く balancing selection の証明ではない。

**両者を同じ生物材料で兼用する必要はない**。Aは技術・表現型軸、Bは機構のアンカー。全球一般則の機構的証明には、成功例だけでなく反証例と系統独立な複数起源を加える必要がある。これらは *pilot* の最低限の試験可能性であり、5仮説を一括して証明する検出力保証ではない。

## H1–H5：新規データがある場合の識別規則

| ID | 一つだけ選ぶ主テスト | 成立するなら | 不成立／未確定の基準 |
|---|---|---|---|
| **H1 phenotype accessibility** | 独立に花器官を反射測色し、事前固定の白―有色軸Wを**粗色ラベル構成保存null**と比較 | 写真外の器官／色校正計測でもWが過剰。少数種の特定の白花率を意味しない | 色齢・露出・萼/苞・二色花の検証で軸が消える、または適格種数不足 |
| **H2 historical sorting** | 同一種の**集団別中立SNP距離**と花色型頻度差を、現在地理距離・気候差と同時比較。利用可能な独立系統地理史があれば追加 | 現在の環境を揃えた別祖先集団ペアで色差が残り、**環境のみの保留予測より中立履歴項が上乗せ** | 中立SNP／独立地史なし＝「history未識別」。距離proxyだけで「分布変遷が優位」としない |
| **H3 environmental sorting** | H2と**同じ集団**で、環境が異なるが中立構造が近い集団ペアを比較。色型×水・光・温度処理の**種子生産/生存適応度交互作用**も測る | 中立系統を調整後の環境効果と適応度交互作用が再現 | 写真のBIO5関連だけ、同一生育環境でのfitnessなし＝選択未証明 |
| **H4 local coexistence maintenance** | 既知の**一つの相互交配集団で個体識別・遺伝型・開花/種子成功**を連続2季以上追跡 | 白と有色が実在する別個体・遺伝型として再出現し、頻度依存又は環境変動に伴う型別fitnessが識別できる | 5km写真ペアや2年の頻度存続だけ＝protected polymorphism未証明。長期安定性にはより長い追跡が必要 |
| **H5 pigment modularity and persistence** | Bと**同じ標識個体**の花弁＋葉でアントシアニンLC-MS/HPLC、候補MYB/F3H/ANSの器官別発現、種子適応度を対応させる | 花弁限定低下（PAL）と葉も低下（WAL）が異なる生理・繁殖反応を示す、複数起源で再現 | PAL/WALが揃わず比較不可、又は白写真だけ＝維持コストの因果機構は未確定。WALの低頻度を欠測値で補完しない |

### 連結パイロットBで同じ個体から採る最小測定

- **全240標識個体（目標）**：標本識別、野生/栽培、色/連続反射スペクトル、花齢、家系可能なら親子関係、中立MIG-seq等のSNP、位置は非公開原票、気温/湿度/土壌、各年の訪花・結実・種子数。
- **組織・発現のサブセット**：2種×2調査個体群×白/有色2型×各6個体＝最大**48個体**の花弁・葉ペア（**96組織**）。WAL型が実在しないときはその比較は**not estimable**であり、勝手な代用品を置かない。MYB/F3H/ANSの配列/オルソログ同定は種別に必要。
- **配置条件**：各種で過去系統または中立遺伝群が異なる地域と現在環境の2水準をできるだけ交差させる。6個体群がすべて同じ中立祖先／同じ気候ならH2/H3の効果は同定できない。実測の独立個体群数に応じて最終判定を縮小する。

### H3/H4の「環境選択 × 繁殖利益」を同じ材料で試す最小操作

H3の現在環境対応とH4の色素型の利益／不利益は、観察された気候相関だけからは分離できない。パイロットBから家系・遺伝型が明らかな繁殖可能な個体を取り、共通圃場で次の小型要因実験を**同時に**行う。

- 要因は **white/chromatic（遺伝的・生化学的に確認）× 水分control/drought × open/supplemental pollen**。型別に各組合せ**独立8個体**を目標とし、2種なら2×2×2×8×2＝**128個体**。実在の別個体／遺伝型がない処理区は欠測とし、存在しなかった型を作らない。集団・母系とブロックを記録して処理を無作為化する。
- 応答は花色反射率、花弁／葉の色素量、植物体サイズ、障害・生存率、受粉・種子数。open − supplemented が訪花不足の指標となり得るが、受粉処理が同一個体の生理・資源状態に与える効果を区別する。
- **最低限の識別**：morph×water（環境ごとにfitness順位が反転するか）と morph×pollen（送粉機会で繁殖利益が変わるか）。独立遺伝背景／家系で反復しなければ色素そのものの費用や利益に帰属できない。
- **測れない主張**：種子差やアントシアニン量の相関だけで「色素合成費用」を金銭的・代謝フラックス的に同定したことにはならない。色素生産コストの因果証明には、遺伝背景を揃えた編集・相補系や安定同位体等の資源流量測定がさらに必要。パイロットでは適応度トレードオフの存在可能性と大きさを推定する。

同じ個体群から実験対象を抽出することでH2/H3の集団史とH4/H5の器官コストに接続できる。各組合せ8個体は**最低限の実行可能性目安**であり、効果の大きさに対する検出力保証ではない。

### 事前解析の数理モデル（パイロットB）

Population \(j\) in species \(i\); color morph frequency \(p_{ij}\) is **genotyped individual frequency**, not an iNaturalist photo frequency.

\[
p_{ij}=\mathrm{logit}^{-1}(\alpha_i+
\beta_E E_{ij}+\beta_H H_{ij}+\beta_G G_{ij}+
u_{i,j})
\]

- \(E\): 水分・花期の温度や訪花構造（可用環境/Mを記録）
- \(H\): 中立SNPに基づく ancestry/集団史指標（独立源の地史情報がなければ history とは呼ばない）
- \(G\): 現在地理距離／個体群構造のモデル化
- 予測比較：**leave-one-population-out**。比較するモデル（geometry only; +E; +neutral ancestry; +both）は、同じ個体群／同じ成功予測行を使用する。自然履歴と現在環境の寄与が共線的なら**分離不可能**と報告。
- H3選択はさらに適応度 \(F\) に対して \(F\sim\mathrm{morph}\times\mathrm{treatment}+\mathrm{lineage}+\mathrm{population}\) を別測定する。「写真の色型頻度→適応度」の代入をしない。
- H4のfrequency-dependent selectionは、稀少型の相対適応度 \(w_{\mathrm{rare}}-w_{\mathrm{common}}\) を複数個体群・複数年で測れなければ「検証不可」。単なるC/S状態分類から推定しない。

**パイロットの成否**はp<.05で決めるだけでなく、予測可能な個体群の有無、同じ環境・異なるancestryペアの有無、少なくとも二つの型を持つ交配集団が現実にあるか、推定対象のCIを併記して判断する。因果性には独立由来情報または操作実験が必要。

## 今ある証拠の「実行済」と「未実行」を分ける

| ID | 2026-10-20実測／公開資料 | 追加の最短ゲート |
|---|---|---|
| H1 | 158種 prospective W = **.517246** vs label-preserving null **.457143** (p=.001)。near-clip制御は **142/158=89.87%**, 事前技術ゲート **INDETERMINATE**。 | 異なる撮影・測色ソース＋独立花器官QC（新規実測必要） |
| H2 | 三つの元コホートで photo phenotypic IBD-like > BIO5 IBE-like（直接neutral SNPは無し）。 | 同一個体群のneutral SNP/歴史障壁（新規または独立アーカイブ必要） |
| H3 | 同種・同暦月・高純度 **327種**のBIO5白−非白は+0.32985℃、そのうち1°小格子間+0.27158℃ (5°null p=.005)、**小格子内+0.05827℃** (小格子内null p=.075)。**約82%は1°セル間の構成差**。これは原因寄与率ではない。 | 中立系譜・地点内処理×形質の適応度（新規実測必要） |
| H4 | 5km同年月・別観察者194種では追加局所類似p=.300。文献厳格C/S 34種とphoto 50km解析は重複5種、**S-onlyの重複0種**。 | 現地の同一交配集団×別個体×2年以上・繁殖適応度 |
| H5 | 先行文献PAL13例/WAL13例（文献選択バイアスあり）。FCP分子登録14種/24比較、独立起源24回ではない。直接FCP写真個体に対応する遺伝型は無し。 | 同一自然集団でpetal・leaf生化学とgenotype/fitnessを同時測定 |

**重要な反証条件**：白色photo軸Wは支持されても、白色型が他色型より強い地理的選択を受けるという仮説は、50kmおよび等面積格子の各成分baseline補正で**非支持**である。この負の結果はH1の反例ではない。

## 提出判断の3段階

1. **現状**：全球における同種内**写真花色の共通軸**＋地域構造＋小さな条件付き気候関連は記述可能。分布史優位／遺伝的FCP維持／組織特異性色素の適応原因は未確定。
2. **A実施後**：独立花器官測色でもWが残れば、写真系への依存を部分的に排除できる。GLOBAL prevalenceの推定にはならない。
3. **B実施後**：歴史と環境を一つの集団で識別する pilot、およびH4/H5の遺伝型・組織・適応度の「繋がった機構」へ初めて進める。2種の結果を全球的原因割合へ外挿しない。

### アーカイブ上の出所

- H1: \`results/polymorphism_h2_third_cohort_prospective_white_axis_20260917/result.json\`; \`results/polymorphism_h2_third_cohort_highlight_validity_20260922/result.json\`.
- H2: \`results/polymorphism_phenotypic_IBD_IBE_posthoc_20261007/result.json\`.
- H3: \`results/fcp_h1h5_minimal_h3_nested_worldclim_verified_20261020.json\`; original-photo CI [#39430364175](https://github.com/zuizui0223/fcp/actions/runs/39430364175).
- H4: \`results/fcp_local_photo_coexistence_month_null_verified_20261018.json\`; \`results/fcp_cs_literature_photo_bridge_20261020/summary.json\` from separate PR #248 source commit 69c0fb481f361fd8764d376a10fbd7035fd73d73.
- H5: \`results/polymorphism_crossspecies_pal_wal_frequency_20260925/result.json\`; \`results/fcp_molecular_expansion_20261019/summary.json\` from the same PR #248 commit.
- Machine/logic verification: \`scripts/analysis/audit_fcp_h1_h5_minimum_20261020.py\`; \`.github/workflows/fcp-h1-h5-minimum-20261020.yml\`. These do not synthesize biological experiments.
