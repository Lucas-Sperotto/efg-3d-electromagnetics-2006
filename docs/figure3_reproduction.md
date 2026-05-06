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

**Decisão desta reprodução:** a figura principal deve usar o caso regular
`15x15x15`. Ele é a malha regular mais refinada, tem a melhor concordância na
janela central e exige muito menos iterações de GMRES. As nuvens `nonuniform`
e `nonuniform_refine` são documentadas como comparações suplementares, não
como figura principal.

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
| GMRES — iterações | 69 |
| GMRES — resíduo relativo | 9,54×10⁻¹⁰ |
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
| `data/output/figure3_region_error_summary.csv` | Comparação por região entre regular e nuvens não uniformes |

Os arquivos de figuras são ignorados pelo `.gitignore`; o CSV de dados e os
arquivos de métricas também. Todos podem ser regenerados executando:

```bash
./build/reproduce_cube_sparse --case plane15
python3 scripts/plot_cube_plane.py --all
```

---

## 4. Métricas principais

### 4.1 Plano completo (0 ≤ y ≤ 10, 0 ≤ z ≤ 10)

| Métrica | Valor |
|---|---|
| Pontos no plano | 10 201 |
| plane_max_abs_error | 8,8854 V |
| plane_mean_abs_error | 0,03976 V |
| plane_relative_error | 7,1030 % |

### 4.2 Interior do plano (0 < y < 10, 0 < z < 10)

| Métrica | Valor |
|---|---|
| Pontos interiores | 9 801 |
| interior_max_abs_error | 4,5832 V |
| interior_mean_abs_error | 0,03399 V |
| interior_relative_error | 4,5519 % |

### 4.3 Centro do cubo (y ≈ 5, z ≈ 5)

| Grandeza | V_num | V_exact | abs_error |
|---|---|---|---|
| Potencial | ≈ 1,638 V | ≈ 1,659 V | ≈ 0,021 V |

O erro absoluto médio no interior (~0,034 V) representa aproximadamente
0,34 % de V₀. Próximo ao centro do plano, o erro fica em ~0,21 % de V₀.

### 4.4 Regiões críticas e comparação com a nuvem não uniforme

Para documentar a região próxima às arestas superiores, foram avaliadas duas
faixas no plano:

- `upper_edge_band_all`: `z >= 9` e (`y <= 1` ou `y >= 9`), incluindo bordas;
- `upper_edge_band_open`: `9 <= z < 10` e (`0 < y <= 1` ou `9 <= y < 10`).

| Caso | Região | Pontos | max abs error | mean abs error | relative error |
|---|---|---:|---:|---:|---:|
| regular 15x15x15 | plano interior | 9801 | 4,5832 V | 0,03399 V | 4,5519 % |
| regular 15x15x15 | upper edge open | 200 | 4,5832 V | 0,66973 V | 18,9098 % |
| regular 15x15x15 | janela central | 3721 | 0,06822 V | 0,01839 V | 0,9304 % |
| nonuniform | plano interior | 9801 | 4,5780 V | 0,04278 V | 5,1413 % |
| nonuniform | upper edge open | 200 | 4,5780 V | 0,80546 V | 20,7600 % |
| nonuniform | janela central | 3721 | 0,09069 V | 0,01942 V | 0,9919 % |
| nonuniform refine | plano interior | 9801 | 4,2374 V | 0,05639 V | 4,4404 % |
| nonuniform refine | upper edge open | 200 | 4,2374 V | 0,62385 V | 17,0825 % |
| nonuniform refine | janela central | 3721 | 0,17985 V | 0,05339 V | 2,8962 % |

A janela central confirma a concordância visual da Figura 3: o erro relativo
fica abaixo de 1 % nos dois casos. A faixa superior concentra a discrepância,
como esperado para uma fronteira descontínua e uma série analítica truncada.

---

## 5. Por que o erro máximo do plano ainda é alto perto da tampa

A política Dirichlet atual torna a fronteira numérica compatível com a série
analítica seno-seno-sinh:

```
V = 0    em x = 0, x = L, y = 0, y = L, z = 0
V = V₀   em z = L e 0 < x < L, 0 < y < L
```

Com essa convenção, as arestas superiores são aterradas e o erro máximo 3D
do caso `15x15x15` cai para `4,642650e-01`, em vez de permanecer próximo de
`V₀`. O plano exportado ainda mostra erro alto em pontos de fronteira próximos
às arestas superiores, por exemplo `y = 0,2`, `z = 10`, onde:

```text
V_num   ≈ 2,975 V
V_exact ≈ 11,860 V
abs_error ≈ 8,885 V
```

Esse ponto combina três dificuldades:

- a série analítica está truncada em 25 termos e exibe sobressinal de Gibbs
  perto da descontinuidade da tampa;
- o EFG/MLS não é interpolante de Kronecker, então pontos de amostragem na
  fronteira entre restrições não precisam reproduzir exatamente `V₀`;
- o gradiente é mais acentuado próximo às arestas superiores.

Por isso, a métrica interior do plano continua sendo a melhor medida para a
comparação visual da Figura 3.

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

- O erro relativo interior do plano obtido aqui (4,55 %) é maior do que o 1,40 %
  reportado pelo artigo para a nuvem não uniforme.
- O erro máximo no interior do plano (4,583 V) é maior do que o 0,15 V
  reportado pelo artigo.
- As iso-linhas próximas às arestas superiores divergem mais do que as
  iso-linhas centrais.

A grade regular é suficiente para a comparação qualitativa da Figura 3, mas
não é suficiente para reproduzir os números quantitativos reportados no artigo.
Na nuvem não uniforme refinada (`base_n=13`, `top_n=15`), o erro 3D interior
fica ligeiramente abaixo do regular `15x15x15` (`2,78 %` contra `2,80 %`) e o
erro relativo interior do plano cai para `4,44 %`. Porém, a janela central
piora para `2,90 %` e o GMRES sobe para `832` iterações. Por isso, a figura
principal permanece sendo o caso regular, enquanto a nuvem refinada é uma
evidência de convergência que precisa de precondicionamento ou melhor
distribuição de pontos.

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

### 7.6 Comparação com os valores reportados no artigo

O artigo reporta erro relativo de aproximadamente `1,40 %` e erro máximo de
aproximadamente `0,15 V`. Esses valores ainda não foram reproduzidos
quantitativamente:

- o caso regular `15x15x15` tem erro relativo 3D interior `2,80 %`;
- o plano regular tem erro relativo interior `4,55 %`;
- o erro máximo 3D amostrado no regular é `0,464 V`;
- o erro máximo interior do plano regular é `4,583 V`, concentrado perto das
  arestas superiores.

As diferenças decorrem de quatro fatores principais: grade regular em vez da
nuvem adaptada do artigo, norma discreta em vez da integral da Eq. (16), série
analítica truncada com Gibbs perto da tampa e ausência de comparação digital
com os níveis exatos da figura publicada.

---

## 8. Conclusão

Esta etapa constitui uma **reprodução qualitativa robusta da Figura 3** do
artigo de Parreira et al. (2006).

O pipeline EFG esparso com GMRES — grade regular 15×15×15, spline cúbica,
domínio de influência cúbico, suporte pelo terceiro vizinho × 1,5, Dirichlet
por multiplicadores de Lagrange — reproduz a distribuição de potencial no plano
x = 5,33 com:

- erro relativo interior do plano de **4,55 %**;
- erro médio interior de **0,034 V** (0,34 % de V₀);
- iso-linhas visualmente concordantes com a solução analítica no interior do
  domínio;
- todos os critérios de segurança satisfeitos (support_lt_4 = 0,
  mls_failures = 0, GMRES convergido).

A concordância qualitativa — forma e posicionamento das linhas de contorno no
interior — é compatível com a Figura 3. A nova política de Dirichlet aterra as
arestas superiores, reduzindo o `max_abs_error` 3D do caso `15x15x15` para
`4,642650e-01`. O erro máximo do plano ainda é dominado por pontos próximos à
tampa e deve ser interpretado junto com as métricas interiores.

**A reprodução quantitativa completa dos resultados do artigo** — erro de
1,40 % e erro máximo de 0,15 V — ainda não foi alcançada. Ela requer uma nuvem
não uniforme mais fiel à distribuição do artigo, com concentração nos cantos
superiores do cubo onde o gradiente é mais acentuado. Uma primeira versão dessa
nuvem foi implementada e está documentada na seção 9 abaixo.

---

## 9. Comparação: grade regular 15×15×15 versus nuvem não uniforme

### 9.1. Configuração da nuvem não uniforme

A função `cube_generate_article_cloud` implementa a seguinte estratégia:

1. **Grade base regular** base_n³: esqueleto uniforme com espaçamento h = L/(base_n−1).
2. **Face superior completa** top_n × top_n: grade em z = L com espaçamento
   h_fine = L/(top_n−1); somente o interior aberto da face recebe V = V₀.
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
| GMRES iterações | 69 | 67 | 391 |
| Erro relativo interior 3D | **2,80 %** | 5,13 % | **2,96 %** |
| max_abs_error 3D | 0,464 V | 0,852 V | 0,834 V |
| Tempo de montagem | 0,67 s | 0,41 s | 0,58 s |
| Tempo de solução | 0,026 s | 0,008 s | 0,150 s |

O refinamento `base_n=13`, `top_n=15` acrescenta:

| Métrica | Nuvem refinada |
| --- | --- |
| Nós totais | 3 089 |
| Restrições | 1 082 |
| support_lt_4 | 0 |
| mls_failures | 0 |
| GMRES iterações | 832 |
| Erro relativo interior 3D | 2,78 % |
| max_abs_error 3D | 0,941 V |
| Tempo de montagem | 0,83 s |
| Tempo de solução | 0,716 s |

### 9.3. Análise comparativa

**Melhora em relação ao 11×11×11 regular:** a nuvem não uniforme (1 974 nós) reduziu o
erro interior de 5,13 % para 2,96 % com apenas 643 nós adicionais em relação ao 11×11×11
regular. Isso confirma que adicionar nós na zona superior melhora a resolução do gradiente
acentuado próximo à tampa z = L.

**Comparação com 15×15×15 regular:** a grade regular com 3 375 nós ainda supera a nuvem não
uniforme por margem pequena (2,80 % vs. 2,96 %). Isso é esperado porque:

- O artigo usa uma nuvem especificamente adaptada à geometria e ao campo de gradiente do
  cubo eletrostático, com concentração de nós nas regiões de gradiente extremo (y → 0 e
  y → L próximos a z = L).
- A nuvem implementada adiciona nós em toda a faixa interior da zona superior, sem
  concentração específica nas arestas onde o gradiente é máximo.
- A métrica de erro usada aqui (L2 discreta em grade uniforme 11×11×11) difere da norma
  contínua do artigo (Eq. 16), tornando a comparação indireta.

**GMRES:** a nuvem não uniforme exigiu 391 iterações contra 69 do 15×15×15 regular. O raio
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
uniformes. O refinamento `base_n=13`, `top_n=15` melhora a métrica 3D interior, mas
expõe um gargalo claro de GMRES; o próximo passo natural é precondicionamento ou uma
distribuição de pontos mais focada nas arestas superiores.
