#import "@preview/physica:0.9.2": *
#import "@preview/cetz:0.3.0"

// simple page setup
#let mainfont = ("BIZ UDPMincho")
#set text(font: mainfont, lang: "ja", size: 9pt)
#set page(numbering: "1")
#set math.equation(numbering: "1.")
#let zeros = $op("Zeros")$
#let poles = $op("Poles")$
#let int = $op("Int")$
#let geigval = $op("geigval")$
#let geigvec = $op("geigvec")$
#let geigpair = $op("geigpair")$
#let neigval = $op("neigval")$
#let neigvec = $op("neigvec")$
#let neigpair = $op("neigpair")$
#let jordancurve = $op("JordanCurve")$
#let dker = $op("dker")$
#let anw = $C (W,CC)$
#let poly = $cal(P)$
#let monicpoly = $cal(M P)$
#let fop = $op("FOP")$
#let fonp = $op("FONP")$
#let mfop = $op("MFOP")$
#let cfop = $op("CFOP")$
#let sfop = $op("SFOP")$
#let hankel = $op("Hankel")$
#let definite = $op("Definite")$
#let gl = $op("GL")$
#let ip(x, y) = $lr(angle.l #x, #y angle.r)$

// regions
#let rd = $RR^d$

== Locating zeros and its multiplicities of an analytic function using numerical integration @kravanja_derivative-free_1999 (@kravanja_locating_1999)

- 行列の添字は1番目が縦方向, 2番目が横方向に行くにつれて増える
- $jordancurve := {im gamma | gamma in C(SS^1,CC) "は単射"}$
- $W$: 単連結, $f in anw (abs(dker f) in NN)$, $gamma in jordancurve (ker f sect gamma = emptyset)$
- $C(W, C)$: $W$上の正則関数の全体, $M(W, CC)$: $W$上の有理関数の全体
- (零点の集合)$ker: anw -> 2^CC$: 零点, $dker -> 2^(CC times NN)$: 零点の多重集合（重複度）
- (零点の個数)$N := abs(dker f), n := abs(ker f)$//$, {(z_i, nu_i)}_{i <= n} := dker f$
- (留数定理)$forall f in M(W, CC). 1/(2 pi i) integral_gamma f(z) dd(z) = sum_z Res(f, z), Res(f, z)$: $f$のローラン展開における$-1$次の係数
- (偏角の原理)$N = 1/(2 pi i) integral_gamma f'/f dd(z)$
- ($f$のみを用いた$N$の推定)$forall c_1, ..., c_G in gamma (arg c_i < arg c_(i+1) and abs(arg_(z = c_k)^c_(k+1) f(z)) < pi). N = 1/(2 pi) sum_(k=1)^G arg f(c_(k+1))/f(c_k)$
  - $because N = 1/(2 pi i) integral_f(gamma) dd(w)/w$
  - ($f$のみを用いた$N$の下限の推定(条件を緩めた場合))$forall c_1, ..., c_G in gamma (arg c_i < arg c_(i+1)). N >= 1/(2 pi) sum_(k=1)^G arg f(c_(k+1))/f(c_k)$

=== Formal orthogonal polynomials
- $poly$: $CC$係数多項式の全体, $poly_t$: $t$次$CC$係数多項式の全体
- $monicpoly$: $CC$係数モニック多項式の全体, $monicpoly_t$: $t$次$CC$係数モニック多項式の全体
- (FOP)$fop: NN times L(poly,CC) -> 2^poly, (t, c) |-> {phi in poly_t | forall i in {0, ..., t-1}. c(z^i phi_t) = 0}$
- (FOPの言い換え,Yule-Walker system)$forall c, t. c_i := c(z^i), sum_(j in {0,...,t}) a_j x^j in fop_(t,c) <==> forall i in {0,..., t-1}. sum_(j in {0,...,t}) a_j c^(i + j) = 0 <==> hankel({c_i}_(i=0)^(t-1)) vec(a_(t-1),dots.v,a_0) = - a_t vec(c_t,dots.v,c_(2t-1))$
- (FOPの行列表記)$forall phi in fop_t. phi prop det mat(c_0,...,c_t;dots.v,dots.down,dots.v;c_(t-1),...,c_(2t-1);1,...,z^t) because "クラメルの公式"$
- (Formal orthogonal polynomialの正規化)
  - (Orthonormality)$forall t, c. fonp_(t,c) := {phi in fop_(t,c) | c(phi^2) = 1}$
  - (Monic)$forall t, c. mfop_(t,c) := {sum_(j in {0,...,t}) a_j x^j in fop_(t,c) | a_t = 1}$
  // - (Constant)$forall t, c. cfop_(t,c) := {sum_(j in {0,...,t}) a_j x^j in fop_(t,c) | a_0 = 1}$
  // - (Sum)$forall t, c. sfop_(t,c) := {sum_(j in {0,...,t}) a_j x^j in fop_(t,c) | sum_(j in {0,...,t}) a_j = 1}$
- (Hankel行列)$forall m in NN. hankel_m: {a_i}_(i=0)^(2m-1) -> CC^(m times m), (i,j) |-> s_(i+j), hankel := union.dot_(m in NN) hankel_m$
- (Definite)$definite_t := {c in L(poly,CC)|det(hankel({c_i}_(i=0)^(t-1))) != 0 <==> hankel({c_i}_(i=0)^(t-1)) in gl(t)}$
- (DefiniteとFOPの存在)$forall c in L(poly,CC). c in definite_t <==> a_k != 0 <==> fop_(t,c) != emptyset$

=== Derivative-dependent method

- 零点の重複度を推定できる, 極とその重複度も実は推定できる
- $c: poly -> CC, phi |-> integral_gamma phi(z) (f'(z))/(f(z)) dd(z) in L(poly,CC)$
- ($c$の意味)$forall phi in poly. c(phi) = sum_((z, nu) in dker f) nu phi(z)$
  - $because poly subset anw therefore poles(phi(z)(f'(z))/f(z)) = poles(1/f(z)) = dker f therefore$留数定理
- (モーメント)この文脈では$c_i$をモーメントという
- (モーメントのHankel行列)$forall m in NN. H_m := hankel({c_i}_(p=0)^(m-1)), H^<_m := hankel({c_i}_(p=1)^(m))$
- (モーメントのHankel行列の零点による表記)$H_m = sum_((z, nu) in dker f) nu_k vec(1,dots.v,z^m) (1,...,z^m)$
- (ヴァンデルモンド行列)$Z_n := diag({z_i}_(i=1)^n), D_n := diag({nu_i}_(i=1)^n), V_n := (i,j) |-> z_i^j in CC^(n times n)$
- (モーメントのHankel行列の対角化)$H_n = V_n D_n V_n^T, H^<_n = V_n D_n Z_n V_n^T$
- (モーメントのHankel行列の$rank$)$forall p in NN_0. rank H_(n + p) = n$
  - $because rank H_(n + p) >= n because H_n "は"H_(n+p)"に含まれる" and H_n in gl(n) because D_n, V_n in gl(n) because z$は互いに異なる
  - $because rank H_(n + p) <= n because "(モーメントのHankel行列の零点による表記)より線型独立な列(行)は最大で"n"個"$
- (FOPの次数の上限)$forall p in NN. fop_(c, n + p) = emptyset$
- ($fop_n$)$fop_n = {z |-> product_(z' in ker f) (z - z')}$

=== Derivative-free method

- 零点の重複度は不明, 極は検出できない


=== Sakurai-Sugiura method @sakurai_projection_2003

- $zeros$: 零点の集合, $geigval(A, B), geigvec(A, B), geigpair(A, B)$: 一般化固有値問題$A v = lambda B v$の一般化固有値の集合,一般化固有ベクトルの集合,2つの組の集合
- $Omega subset CC$: 単連結, $F(z): Omega -> CC^(n times n) (not forall z in Omega. det F(z) = 0 and "正則")$
- (非線形固有値問題)$neigval(F), neigvec(F), neigpair(F)$: $F(lambda) x = 0$を満たす$lambda$の集合, $x$の集合, 2つの組の集合
- (単因子標準形)$exists P, Q: Omega -> CC^(n times n) (exists p, q in CC. forall z in Omega. det P(z)= p and det Q(z) = q). forall j in {1, ..., n}. exists d_j (z) (forall j in {2, ..., n}. (d_j (z))/(d_(j-1) (z)) "は正則"). exists D: Omega -> CC^(n times n), z |-> diag(d_1 (z), ..., d_n (z)). forall z in Omega. P(z)F(z)Q(z)=D(z)$
- ${lambda_1, ..., lambda_s} := zeros(d_n)$, $p_j$: $P$の行ベクトルの転置, $q_j$: $Q$の列ベクトル
- $forall j in {1, ..., n}. exists h_j: Omega -> CC (forall z in Omega. h_j (z) != 0, h_j "は解析的"). forall i in {1, ..., s}. exists alpha_(j i) in NN (alpha_(j-1 i) <= alpha_(j i)). d_j (z) = h_j (z) product_(i in {1, ..., s}) (z - lambda_i)^(alpha_(j i))$
- ($d$の零点と固有値・固有ベクトル)$forall j in {1, ..., n}. forall lambda_i in zeros(d_j). (lambda_i, q_j (lambda_i)) in neigpair(F)$
  - $because forall j in {1, ..., n}. forall lambda_i in zeros(d_j). F(lambda_i) q_j (lambda_i) = P(lambda_i)^(-1) D(lambda_i) Q(lambda_i)^(-1) (Q(lambda_i) e_j) = P(lambda_i)^(-1) (D(lambda_i) e_j) = P(lambda_i)^(-1) d_j (lambda_i) e_j = 0 because D$は対角行列
  - $because d_j (lambda_i) = 0 (because lambda_i in zeros(d_j))$
- $forall u, v in CC^n. chi_j (z) := u^T q_j (z) p_j (z)^T v. u^T F(z)^(-1) v = sum_(j in {1, ..., n}) (chi_j (z))/(d_j (z))$
- (モーメント)$forall k in NN_0. mu_k := 1/(2 pi i) integral_Gamma z^k u^T F(z)^(-1) v dd(z)$
  - $mu_k = 1/(2 pi i) integral_Gamma z^k f(z) dd(z) = 1/(2 pi i) integral_Gamma z^k sum_(j = 1)^n (chi_j (z))/(d_j (z)) dd(z) = sum_(l = 1)^m (chi_n (lambda_l))/(d'_n (lambda_l)) lambda_l^k = sum_(l = 1)^m nu_l lambda_l^k, nu_l := (chi_n (lambda_l))/(d'_n (lambda_l))$
- (モーメント2)$forall k in NN_0. s_k := 1/(2 pi i) integral_Gamma z^k F(z)^(-1) v dd(z)$ // NN_0でよい？
  - $s_k = sum_(l=1)^m sigma_l lambda_l^k, sigma_l := (q_n (z) p_n (z)^T v)/(d'_n (z))$
- $forall Gamma (int Gamma subset Omega).$, $lambda_1, ..., lambda_m$: $Gamma$内の異なる縮退しない固有値
- (ハンケル行列)$H_m := (mu_(i+j-2))_(i,j in {1, ..., m}), H^<_m := (mu_(i+j - 1))_(i,j in {1, ..., m})$
- (Sakurai-Sugiura method, 一般化固有値)$forall l in {1, ..., m}. chi_n (lambda_l) != 0 ==> {lambda_1, ..., lambda_m} = geigval(H^<_m, H_m)$
  - $therefore V_m := mat(1,1,...,1;lambda_1,lambda_2,...,lambda_m;dots.v,dots.v,space,dots.v;lambda_1^(m-1),lambda_2^(m-1),...,lambda_m^(m-1)), D_m := diag(nu_1, ..., nu_m), Lambda_m := diag(lambda_1, ..., lambda_m), H^<_m - lambda_l H_m = V_m D_m (Lambda_m - lambda_l I) V_m^T and forall l in {1, ..., m}. nu_l != 0 (because chi_n (lambda_l) != 0)$
- $S := (s_i)_(i in {0, ..., m - 1})$
- $forall (lambda_l, w_l) in geigpair(H^<_m, H_m). exists c_l in CC without {0}. q_n (lambda_l) = c_l S w_l$
  - $because forall (lambda_l, w_l) in geigpair(H^<_m, H_m). (H^<_m - lambda_l H_m) w_l = 0 therefore V_m D_m (Lambda_m - lambda_l I) V_m^T w_l therefore (Lambda_m - lambda_l I) V_m^T w_l because V_m, D_m$は正則
  - $therefore exists alpha_l in CC without {0}. V_m^T w_l = alpha_l e_l because (Lambda_m - lambda_l I) = diag(lambda_1 - lambda_l, ..., lambda_(l-1) - lambda_l, 0, lambda_(l+1) - lambda_l, ... lambda_m - lambda_l) and lambda_1, ..., lambda_m$は相異なる
  - $therefore S = (s_0, ..., s_(m-1)) = (beta_1 q_n (lambda_1), ..., beta_m q_n (lambda_m)) V_m^T, beta_l := (p_n (lambda_l)^T v)/(d'_n (lambda_l))$
  - $therefore q_n (lambda_l) = (S V_m e_l) beta_l = (S w_l)/(alpha_l beta_l)$
- (Sakurai-Sugiura method, 一般化固有ベクトル)$forall (lambda_j, w_j) in geigpair(H^<_m, H_m). (lambda_j, q_n (lambda_j)) in neigpair(F)$

=== Circle case

- $Gamma := {gamma + rho e^(i theta) | theta in [0, 2 pi)}$
- $mu_k approx 1/N sum_(j = 0)^(N - 1) ((omega_j - gamma)/rho)^(k + 1) f(omega_j), omega_j := gamma + rho e^(2 pi i (j + 1/2)/N)$

=== Block Sakurai-Sugiura method

- $f$をバッチで計算する



#bibliography("main.bib")
