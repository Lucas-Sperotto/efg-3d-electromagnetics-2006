# Reprodução qualitativa da Figura 3

**Referência:** Parreira et al., "The Element-Free Galerkin Method in
Three-Dimensional Electromagnetic Problems", *IEEE Transactions on Magnetics*,
vol. 42, n. 4, pp. 711–714, 2006.

---

## 1. Objetivo da figura

A Figura 3 do artigo exibe as linhas de contorno do potencial elétrico no plano
*yz* em x = 5,33, obtidas tanto pelo método EFG quanto pela solução analítica.
O plano x = 5,33 foi escolhido por ser aproximadamente o plano médio do cubo
L = 10, suficientemente afastado das paredes laterais para evidenciar o
gradiente vertical do potencial entre o fundo aterrado (z = 0) e a tampa em
V₀ = 10 V (z = L).

A figura demonstra que as linhas de contorno do EFG coincidem visualmente com
as linhas analíticas no interior do domínio, mesmo com os fortes gradientes
próximos aos cantos superiores da caixa.

---

## 2. Configuração numérica

| Parâmetro | Valor |
|---|---|
| Problema | Cubo eletrostático L = 10, V₀ = 10 V |
| Distribuição de nós | Grade regular 15×15×15 (3 375 nós) |
| Células de integração | 15×15×15 |
| Quadratura | Gauss-Legendre 2×2×2 (8 pontos/célula) |
| Função peso | Spline cúbica |
| Domínio de influência | Cúbico |
| Raio de suporte | 3° vizinho mais próximo × 1,5 |
| Condições de Dirichlet | Multiplicadores de Lagrange |
| Armazenamento | Esparso COO → CSR |
| Solver | GMRES reiniciado (restart = 300, tol = 10⁻⁹) |
| Nós com suporte ativo < 4 | 0 |
| Falhas MLS (solução) | 0 |
| Falhas MLS (plano) | 0 |
| GMRES — iterações | 70 |
| GMRES — resíduo relativo | 9,03×10⁻¹⁰ |
| Plano exportado | x = 5,33 |
| Grade de amostragem no plano | 101 × 101 pontos (y × z) |

A solução de referência é a série analítica truncada com 25 termos por dimensão
(n, m = 1, 3, 5, …, 49), conforme as Eqs. (13)–(15) do artigo.

---

## 3. Arquivos gerados

| Arquivo | Conteúdo |
|---|---|
| `data/output/cube_plane_x_5_33_refine15.csv` | Dados do plano: x, y, z, V_num, V_exact, abs_error (10 201 pontos) |
| `data/output/figures/cube_plane_x_5_33_overlay_isolines.png` | Sobreposição de iso-linhas EFG (azul) e analítica (vermelho tracejado) |
| `data/output/figures/cube_plane_x_5_33_article_style.png` | Figura no estilo do artigo: iso-linhas V = 1, 2, ..., 10 V |
| `data/output/figures/cube_plane_x_5_33_V_num_isolines.png` | Iso-linhas do potencial numérico isolado |
| `data/output/figures/cube_plane_x_5_33_V_exact_isolines.png` | Iso-linhas do potencial analítico isolado |
| `data/output/figures/cube_plane_x_5_33_V_num_contour.png` | Mapa de contorno preenchido — V_num |
| `data/output/figures/cube_plane_x_5_33_V_exact_contour.png` | Mapa de contorno preenchido — V_exact |
| `data/output/figures/cube_plane_x_5_33_abs_error_contour.png` | Mapa de contorno preenchido — |V_num − V_exact| |
| `data/output/figures/cube_plane_x_5_33_comparison.png` | Painel comparativo com os três mapas preenchidos |
| `data/output/figures/cube_plane_x_5_33_metrics.txt` | Métricas completas do plano |
| `data/output/figures/cube_plane_x_5_33_metrics.csv` | As mesmas métricas em formato CSV |

Os arquivos de figuras são ignorados pelo `.gitignore`; o CSV de dados e os
arquivos de métricas também. Todos podem ser regenerados executando:

```bash
./build/reproduce_cube_sparse --case plane15
python3 scripts/plot_cube_plane.py
```

---

## 4. Métricas principais

### 4.1 Plano completo (0 ≤ y ≤ 10, 0 ≤ z ≤ 10)

| Métrica | Valor |
|---|---|
| Pontos no plano | 10 201 |
| plane_max_abs_error | ≈ 10,00 V |
| plane_mean_abs_error | 0,04249 V |
| plane_relative_error | 8,1894 % |

### 4.2 Interior do plano (0 < y < 10, 0 < z < 10)

| Métrica | Valor |
|---|---|
| Pontos interiores | 9 801 |
| interior_max_abs_error | 4,366 V |
| interior_mean_abs_error | 0,03368 V |
| interior_relative_error | 4,5826 % |

### 4.3 Centro do cubo (y ≈ 5, z ≈ 5)

| Grandeza | V_num | V_exact | abs_error |
|---|---|---|---|
| Potencial | ≈ 1,672 V | ≈ 1,660 V | ≈ 0,012 V |

O erro absoluto médio no interior (~0,034 V) representa aproximadamente
0,34 % de V₀. Próximo ao centro do plano, o erro cai para ~0,12 % de V₀.

---

## 5. Por que o erro máximo é ≈ V₀ nas bordas e cantos superiores

O problema do cubo impõe duas condições de Dirichlet que se contradizem nas
arestas superiores:

```
V = 0    em y = 0, y = L, z = 0  (paredes laterais e fundo)
V = V₀   em z = L               (tampa superior)
```

Nos pontos (y = 0, z = L) e (y = L, z = L), as duas condições coexistem sem
solução única. Os dois métodos resolvem essa ambiguidade de maneiras
matematicamente diferentes:

**EFG com multiplicadores de Lagrange:** a condição da tampa superior
(z = L, V = V₀) é aplicada via multiplicadores a todos os nós da face z = L,
incluindo aqueles que pertencem simultaneamente às paredes laterais. Portanto:

```
V_num(y = 0, z = L) ≈ V₀ = 10 V
```

**Série analítica truncada:** a base é um produto de senos,
sin(nπy/L) × sin(mπx/L), que se anula sempre que y = 0 ou y = L,
independentemente de z. Portanto:

```
V_exact(y = 0, z = L) = 0 V   (sin(0) = 0 em todos os termos)
```

A diferença é inevitável:

```
abs_error(y = 0, z = L) = |V₀ − 0| = V₀ = 10 V
```

Esse erro não indica falha no solver. Ele é o reflexo de uma descontinuidade
de condição de contorno que os dois métodos tratam por convenção diferente.

**Fenômeno de Gibbs na série truncada:** próximo à aresta (y ≈ 0,2, z = L),
a série analítica com 25 termos apresenta sobressinal, retornando
V_exact ≈ 11,86 V > V₀ = 10 V. Nesse ponto específico, o EFG — que impõe
exatamente V = V₀ via Lagrange — é numericamente mais correto do que a própria
série de referência. Esse sobressinal afeta apenas os pontos na fronteira
z = L e não contamina as métricas do interior.

**Gradiente acentuado próximo às arestas superiores:** no interior do plano,
o erro máximo (4,366 V) ocorre em y ≈ 0,1, z ≈ 9,8 — o ponto mais próximo
da aresta superior ainda dentro do filtro interior. Essa região tem gradiente
de potencial extremamente acentuado, difícil de resolver com a densidade de
nós atual. O artigo original usa nuvem não uniforme com concentração de nós
justamente nessa região, o que reduz esse erro significativamente.

---

## 6. Concordância no interior

A figura de sobreposição de iso-linhas
(`cube_plane_x_5_33_overlay_isolines.png`) exibe as linhas de nível de
V_num (azul) e V_exact (vermelho tracejado) no mesmo plano, com os níveis
V = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 V.

**Observações sobre a sobreposição:**

- Para V ≤ 6 V (linhas no terço inferior e médio do domínio), as iso-linhas
  EFG e analítica são praticamente coincidentes.
- Para V ≥ 7 V (linhas que se curvam em direção à tampa superior), a
  concordância permanece boa no centro do domínio (y ∈ [2, 8]) e se deteriora
  progressivamente à medida que y → 0 ou y → L.
- O erro médio interior de 0,034 V (0,34 % de V₀) é consistente com o desvio
  visual de poucos pixels entre as linhas na região central.

Essa concordância qualitativa — linhas de contorno sobrepostas no interior —
é o equivalente visual do que a Figura 3 do artigo demonstra.

---

## 7. Limitações

### 7.1 Grade regular versus nuvem não uniforme

O artigo usa uma distribuição de nós não uniforme com concentração nas regiões
de gradiente acentuado (cantos superiores do cubo). Essa estratégia reduz o
erro especificamente onde a solução é mais difícil de resolver. A grade regular
15×15×15 usada aqui tem espaçamento uniforme h = L/14 ≈ 0,714 em todas as
direções e não adapta a resolução ao gradiente local. Como consequência direta:

- O erro relativo interior obtido aqui (4,58 %) é maior do que o 1,40 %
  reportado pelo artigo para a nuvem não uniforme.
- O erro máximo no interior do plano (4,366 V) é maior do que o 0,15 V
  reportado pelo artigo.
- As iso-linhas próximas às arestas superiores divergem mais do que as
  iso-linhas centrais.

A grade regular é suficiente para a comparação qualitativa da Figura 3, mas
não é suficiente para reproduzir os números quantitativos reportados no artigo.

### 7.2 Série analítica truncada como referência

A solução de referência V_exact é calculada por série de Fourier com 25 termos
por dimensão — não é a solução exata infinita. Perto das descontinuidades de
condição de contorno (z = L com y → 0 ou y → L), a série truncada exibe
sobressinal de Gibbs (V_exact > V₀ em alguns pontos da borda z = L). A métrica
abs_error nessas sub-regiões reflete parcialmente o erro da série de referência,
não exclusivamente o erro do EFG. As métricas do interior (0 < y < L,
0 < z < L) são confiáveis porque, nessa região, o máximo da série truncada
fica abaixo de V₀.

### 7.3 Níveis de contorno aproximados

Os níveis de iso-linha usados nas figuras (V = 1, 2, …, 10 V) são uma
estimativa baseada na análise visual do intervalo de potencial. O artigo não
especifica explicitamente os níveis de contorno usados na Figura 3. A
concordância visual depende dos níveis escolhidos; se o artigo usa níveis
diferentes, a comparação deve ser atualizada.

### 7.4 Métrica de erro diferente

O artigo mede o erro pela integral contínua da Eq. (16):

$$
e = \left[
\frac{\int_{\Omega}(V_\text{EFG} - V_\text{exato})^2\,d\Omega}
{\int_{\Omega}V_\text{exato}^2\,d\Omega}
\right]^{1/2}
$$

Este projeto usa a norma discreta análoga calculada sobre uma grade de
amostragem uniforme. Para distribuições de nós uniformes e grades de
amostragem densas, as duas métricas convergem; ainda assim, os valores
numéricos não são diretamente comparáveis sem verificação adicional.

### 7.5 Comparação visual, não pixel a pixel

Esta é uma comparação qualitativa: verificação de que a topologia e o
espaçamento das linhas de contorno no interior do plano são compatíveis com os
do artigo. Não foi realizada extração de iso-linhas da Figura 3 original para
comparação digital.

---

## 8. Conclusão

Esta etapa constitui uma **reprodução qualitativa robusta da Figura 3** do
artigo de Parreira et al. (2006).

O pipeline EFG esparso com GMRES — grade regular 15×15×15, spline cúbica,
domínio de influência cúbico, suporte pelo terceiro vizinho × 1,5, Dirichlet
por multiplicadores de Lagrange — reproduz a distribuição de potencial no plano
x = 5,33 com:

- erro relativo interior de **4,58 %**;
- erro médio interior de **0,034 V** (0,34 % de V₀);
- iso-linhas visualmente concordantes com a solução analítica no interior do
  domínio;
- todos os critérios de segurança satisfeitos (support_lt_4 = 0,
  mls_failures = 0, GMRES convergido).

A concordância qualitativa — forma e posicionamento das linhas de contorno no
interior — é compatível com a Figura 3. O erro máximo absoluto (~V₀ nas arestas
superiores) é estrutural e matematicamente explicado pela descontinuidade de
condição de contorno; não indica falha no método.

**A reprodução quantitativa completa dos resultados do artigo** — erro de
1,40 % e erro máximo de 0,15 V — requer a implementação da distribuição não
uniforme de nós do artigo, com concentração nos cantos superiores do cubo onde
o gradiente é mais acentuado. Uma primeira versão dessa nuvem foi implementada
e está documentada na seção 9 abaixo.

---

## 9. Comparação: grade regular 15×15×15 versus nuvem não uniforme

### 9.1. Configuração da nuvem não uniforme

A função `cube_generate_article_cloud` implementa a seguinte estratégia:

1. **Grade base regular** base_n³: esqueleto uniforme com espaçamento h = L/(base_n−1).
2. **Face superior completa** top_n × top_n: grade em z = L com espaçamento
   h_fine = L/(top_n−1), todos os nós tornam-se Dirichlet V = V₀.
3. **Fatias interiores** na zona z > L × (1 − z_frac): n_extra_slices fatias
   com espaçamento uniforme em z, adicionando apenas os nós interiores em xy
   (ix = 1..top_n−2, iy = 1..top_n−2) para evitar nós extras de superfície lateral.

Parâmetros executados: base_n = 11, top_n = 13, n_extra_slices = 4, z_frac = 0,30.

### 9.2. Comparação de métricas

| Métrica | 15×15×15 regular | 11×11×11 regular | Nuvem não uniforme |
| --- | --- | --- | --- |
| Nós totais | 3 375 | 1 331 | 1 974 |
| Restrições | 1 178 (34,9 %) | 602 (45,2 %) | 762 (38,6 %) |
| support_lt_4 | 0 | 0 | 0 |
| mls_failures | 0 | 0 | 0 |
| GMRES convergiu | YES | YES | YES |
| GMRES iterações | 70 | — | 428 |
| Erro relativo interior | **2,81 %** | 5,55 % | **4,63 %** |
| max_abs_error | ~V₀ | ~V₀ | ~V₀ |
| Tempo de montagem | 0,64 s | — | 0,52 s |
| Tempo de solução | 0,025 s | — | 0,164 s |

### 9.3. Análise comparativa

**Melhora em relação ao 11×11×11 regular:** a nuvem não uniforme (1 974 nós) reduziu o
erro interior de 5,55 % para 4,63 % com apenas 643 nós adicionais em relação ao 11×11×11
regular. Isso confirma que adicionar nós na zona superior melhora a resolução do gradiente
acentuado próximo à tampa z = L.

**Comparação com 15×15×15 regular:** a grade regular com 3 375 nós ainda supera a nuvem não
uniforme (2,81 % vs. 4,63 %). Isso é esperado porque:

- O artigo usa uma nuvem especificamente adaptada à geometria e ao campo de gradiente do
  cubo eletrostático, com concentração de nós nas regiões de gradiente extremo (y → 0 e
  y → L próximos a z = L).
- A nuvem implementada adiciona nós em toda a faixa interior da zona superior, sem
  concentração específica nas arestas onde o gradiente é máximo.
- A métrica de erro usada aqui (L2 discreta em grade uniforme 11×11×11) difere da norma
  contínua do artigo (Eq. 16), tornando a comparação indireta.

**GMRES:** a nuvem não uniforme exigiu 428 iterações contra 70 do 15×15×15 regular. O raio
de suporte máximo passou de 27 para 53 nós ativos, indicando maior variabilidade na
conectividade local — o sistema fica menos uniformemente condicionado sem precondicionador.

### 9.4. O que falta para reproduzir os 1,40 % do artigo

1. **Concentração de nós nas arestas superiores:** o artigo concentra nós especificamente
   em y ≈ 0 e y ≈ L para z próximo de L, onde o gradiente é máximo. A nuvem atual adiciona
   nós no interior em xy (y distante das paredes), que não resolve o gradiente de borda.
2. **Possivelmente mais nós totais:** o artigo não especifica explicitamente o número de
   nós da nuvem não uniforme.
3. **Norma de erro diferente:** 1,40 % é calculado pela integral contínua da Eq. (16);
   nossa norma discreta não é diretamente comparável.

A nuvem atual demonstra melhora qualitativa e valida o pipeline para distribuições não
uniformes. O refinamento para atingir os valores quantitativos do artigo é um próximo passo
identificado no `TODO.md` §11.3.
