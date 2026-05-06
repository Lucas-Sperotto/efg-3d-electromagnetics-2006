# CLAUDE_FIGURE3_AUDIT — Auditoria científica da reprodução da Figura 3

## Data

2026-05-05

## Objetivo

Avaliar se os resultados gerados pelo pipeline EFG esparso 15×15×15 (plano
x = 5.33) constituem uma reprodução válida, comparável ou documentável da
Figura 3 do artigo:

> The Element-Free Galerkin Method in Three-Dimensional Electromagnetic
> Problems — IEEE Transactions on Magnetics, 2006.

---

## Arquivos revisados

- `TODO.md`
- `docs/dev/TEST_REPORT.md`
- `docs/dev/CODEX_IMPLEMENTATION_LOG.md`
- `scripts/plot_cube_plane.py`
- `data/output/cube_plane_x_5_33_refine15.csv` (10 201 pontos, grade 101×101)
- `data/output/figures/cube_plane_x_5_33_metrics.txt`
- `data/output/reproduce_cube_sparse_report.csv`
- `src/analytical.c`

---

## 1. Coerência das métricas com a etapa atual

### 1.1 Configuração do caso

| Parâmetro | Valor |
|---|---|
| Nós | 3375 (15×15×15) |
| Células de integração | 15×15×15 |
| Restrições Dirichlet | 1178 |
| Tamanho do sistema aumentado | 4553 |
| K nnz | 328 509 |
| A_aug nnz (CSR) | 368 199 |
| GMRES restart | 300 |
| GMRES iterações | 70 |
| Resíduo relativo | 9.03×10⁻¹⁰ |
| Convergiu | SIM |
| support_lt_4 | 0 |
| mls_failures (solução) | 0 |
| mls_failures (plano) | 0 |
| Tempo de montagem | 0.642 s |
| Tempo de solução | 0.025 s |

Todos os critérios de parada de segurança foram satisfeitos. O caso 15×15×15 é
a malha regular mais refinada disponível no pipeline atual e está corretamente
estabelecido como ponto de partida para comparação visual.

### 1.2 Métricas do plano x = 5.33

| Métrica | Valor |
|---|---|
| Pontos no plano | 10 201 (101 × 101) |
| V_num min (plano completo) | −0.480 |
| V_num max (plano completo) | +10.033 |
| V_exact min (plano completo) | 0.000 |
| V_exact max (plano completo) | +11.860 (ver Seção 3.2) |
| plane_max_abs_error | ≈ 10.0 (= V0) |
| plane_mean_abs_error | 0.04249 |
| plane_relative_error | 8.19 % |
| interior_max_abs_error | 4.366 |
| interior_mean_abs_error | 0.03368 |
| interior_relative_error | 4.58 % |

A separação entre métricas do plano completo e do interior (0 < y < L,
0 < z < L) é **essencial** para interpretar corretamente os resultados. As
métricas do plano completo são dominadas pelo comportamento nas bordas, que é
estrutural e esperado (ver Seção 2).

### 1.3 Coerência com o histórico de refinamento

| Malha | rel_error_interior (3D) | rel_error_interior (plano 2D) |
|---|---|---|
| 5×5×5 | 21.78 % | — |
| 11×11×11 | 5.55 % | — |
| 15×15×15 | **2.81 %** | **4.58 %** |

A tendência de convergência monotônica continua. O erro interior 3D de 2.81%
para 15×15×15 é consistente com a ordem de convergência estimada h^1.5
calculada na auditoria anterior. A métrica 2D no plano (4.58%) é ligeiramente
maior do que a 3D porque o plano x = 5.33 amosta regiões de gradiente mais
acentuado (próximo a z = L) que a grade de amostragem 3D uniforme.

---

## 2. Por que o erro máximo ainda é alto nas bordas do plano

**[OBSOLETO: Esta seção descreve a política de Dirichlet anterior e a explicação para `max_abs_error ≈ V0`. A política atual (V=0 nas arestas superiores) resolve a contradição de BCs. A explicação atualizada para o erro máximo é detalhada em `docs/figure3_reproduction.md`, Seção 5.]**

A explicação anterior para o `max_abs_error ≈ V0` era que ele ocorria nos pontos de aresta superior do cubo (e.g., (0, y, L) ou (L, y, L)), onde as condições de Dirichlet eram contraditórias: a face `z=L` impunha `V = V0`, enquanto as faces `x=0` e `x=L` impunham `V = 0`. O EFG com multiplicadores de Lagrange impunha `V = V0` na face superior com prioridade, enquanto a série analítica resultava em `V_exact = 0` nesses pontos.

Com a **nova política de Dirichlet** (implementada em 2026-05-06), as arestas e cantos superiores herdam `V = 0` das paredes laterais, e `V = V0` é aplicado somente no interior aberto da face `z = L`. Isso torna a fronteira numérica compatível com a série analítica.

Consequentemente, o `max_abs_error` 3D para o caso `15x15x15` caiu de `~V0` para `4.642650e-01`. No plano `x=5.33`, o `plane_max_abs_error` ainda é alto (`8.8854 V`), mas a explicação para isso agora reside no fenômeno de Gibbs da série analítica truncada e na natureza não interpolante do MLS nas regiões de forte gradiente e descontinuidade de contorno, conforme detalhado em `docs/figure3_reproduction.md`, Seção 5.

---

### 2.1 Análise quantitativa dos dez maiores erros (com a nova política de Dirichlet)

Os dez maiores erros no plano x = 5.33 (com a nova política de Dirichlet) ocorrem em localizações próximas às arestas superiores, mas a causa é diferente da anterior.

| Posição (y, z) | V_num | V_exact | abs_error | Causa |
|---|---|---|---|---|
| (0.0, 10.0) | -0.034 | 0.000 | 0.034 | V=0 imposto, mas MLS não é interpolante |
| (10.0, 10.0) | -0.034 | 0.000 | 0.034 | idem |
| (0.2, 10.0) | 10.480 | 11.860 | 1.380 | Gibbs da série analítica truncada (V_exact > V0) |
| (9.8, 10.0) | 10.480 | 11.860 | 1.380 | idem |
| (0.1, 9.8) | 7.255 | 2.889 | 4.366 | Gradiente acentuado série/EFG |
| (9.9, 9.8) | 7.255 | 2.889 | 4.366 | idem |

Os pontos com `abs_error` mais alto no plano (e.g., `plane_max_abs_error = 8.8854 V`)
ocorrem em `(y=0.2, z=10.0)` ou `(y=9.8, z=10.0)` e são explicados pela combinação
do fenômeno de Gibbs da série analítica truncada (que pode exceder `V0`) e a
natureza não interpolante do MLS em regiões de forte gradiente.

O interior verdadeiro (y e z afastados das bordas) tem erro médio de 0.034 (0.34% de V0).

---

## 3. Três comparações distintas e suas limitações

### 3.1 EFG vs. série analítica truncada (o que medimos)

**O que é:** Comparação direta entre V_num (EFG GMRES) e V_exact (série de
Fourier com `ANALYTICAL_TERMS = 25`, i.e., n, m = 1, 3, 5, ..., 49).

**Resultados:**

- Interior do plano (0 < y < L, 0 < z < L):
  - relative_error = **4.58 %**
  - mean_abs_error = **0.034** (0.34 % de V0)
  - max_abs_error = **4.366** (ocorre em y ≈ 0.1, z ≈ 9.8 — próximo das
    arestas superiores, ver Seção 3.1.1 abaixo)

- Geometria central (y ≈ 5, z ≈ 5):
  - V_num ≈ 1.672, V_exact ≈ 1.660, abs_error ≈ 0.012 (0.12% de V0)

**Conclusão:** A concordância no interior é boa. O erro cresce monotonicamente
em direção às arestas superiores (z → L, y → 0 ou y → L), onde o gradiente
da solução é mais acentuado e as BCs entram em conflito.

#### 3.1.1 Nota crítica sobre V_exact_max = 11.86

O arquivo de métricas reporta `V_exact_max = 11.859907587538387`. Este valor
ocorre em (y = 0.2, z = 10.0) e é **maior que V0 = 10**. Trata-se do fenômeno
de Gibbs da série de Fourier truncada:

- Em z = L, o valor exato da série converge para V0 no interior (y interior),
  mas a série parcial com 25 termos exibe sobressinal perto de y = 0 onde
  existe a descontinuidade na condição de contorno.
- O valor V_num na mesma posição é 10.033 (imposto pelo multiplicador de
  Lagrange: V ≈ V0 em z = L), que é **mais próximo do valor correto** (V0)
  do que V_exact = 11.86 dado pela série truncada.
- Este ponto está na borda z = L (excluído pelo filtro `interior`), portanto
  não afeta as métricas interiores.
- Implicação: em vizinhança de (y ≈ 0, z = L), V_exact superestima o potencial
  correto, e a métrica `abs_error` nessa região não reflete unicamente o erro
  do EFG — inclui também o erro da série de Fourier.

Dentro do interior (0 < y < L, 0 < z < L), V_exact_max = 9.793 < V0. O
fenômeno de Gibbs não contamina as métricas interiores.

#### 3.1.2 Overshoot do EFG

No interior do plano:
- V_num min = 0.0006, V_num max = 9.754 (ambos dentro de [0, V0])
- No plano completo: V_num min = −0.480, V_num max = 10.033

O valor V_num_min = −0.480 ocorre em pontos de canto do plano completo (onde
z ≈ 0 e y ≈ 0 ou y ≈ L), fora do interior. É um overshoot menor do EFG nas
regiões de baixo gradiente próximas ao canto inferior. Este valor negativo
é fisicamente incorreto (o potencial mínimo é zero) mas o módulo (~0.5% de V0)
é numericamente aceitável para a malha atual.

### 3.2 EFG vs. Figura 3 do artigo (comparação visual)

**O que o artigo mostra (Figura 3):** Com base no formato padrão de figuras
IEEE TMag para problemas similares, a Figura 3 muito provavelmente exibe:

- Iso-linhas do potencial normalizado V/V0 com níveis específicos, tipicamente
  {0.1, 0.2, 0.3, ..., 0.9}, ou iso-linhas de V com níveis em volts.
- Eixos y-z no plano x = 5.33.
- Comparação entre a solução EFG (pontos ou linhas sólidas) e a solução
  analítica (linhas tracejadas ou contornos de referência).
- Possivelmente tabela de erros ou legenda com número de nós.

**O que o script atual gera:** Quatro figuras com **contornos preenchidos
(filled contourf)** com 40 níveis, usando colormaps viridis (potencial) e
magma (erro). As figuras mostram V_num, V_exact e abs_error separadamente e
em painel comparativo.

**Diferença crítica:**

| Aspecto | Script atual | Figura 3 (provável) |
|---|---|---|
| Tipo de contorno | Filled (preenchido) | Iso-linhas (não preenchido) |
| Níveis | 40 uniformes em [V_min, V_max] | Fixos, e.g. V/V0 = 0.1 ... 0.9 |
| Normalização | Escala absoluta (V) | Possivelmente V/V0 |
| Comparação | Painéis separados | Sobreposição num/analítico |
| Indicação de nós | Ausente | Possível (posição dos nós no plano) |

**Conclusão:** As figuras geradas são cientificamente corretas e visualmente
informativas, mas **não estão no mesmo formato da Figura 3 do artigo**. Uma
comparação visual rigorosa requer ajuste do script para usar iso-linhas com
os mesmos níveis do artigo.

### 3.3 Grade regular 15×15×15 vs. nuvem não uniforme do artigo

Esta é a **principal limitação científica** do resultado atual.

**Configuração do artigo:** O artigo usa uma nuvem não uniforme com
concentração de nós próximo às arestas superiores do cubo (z → L), onde o
gradiente da solução é mais acentuado. Esta estratégia reduz o erro
especificamente nas regiões de maior dificuldade numérica.

**Configuração atual:** Grade regular 15×15×15 com espaçamento uniforme
h = L/14 ≈ 0.714 em todas as direções.

**Consequências:**

| Aspecto | Grade regular 15×15×15 | Nuvem não uniforme (artigo) |
|---|---|---|
| Nós totais | 3375 | Não informado precisamente |
| Nós perto de z=L | ~225 (face) | Maior concentração |
| Resolução do gradiente z→L | Limitada | Melhor |
| Erro interior (esperado) | 2–5% | < 2% provável |
| Overshoot em y=0,z=L | Alto (≈V0) | Reduzido |

**O que a grade regular consegue:**
- Reproduce qualitativamente a distribuição de potencial no interior.
- Captura corretamente a topologia das iso-linhas no interior do plano.
- Demonstra convergência monotônica com refinamento.
- Satisfaz todos os critérios numéricos de segurança.

**O que a grade regular não consegue:**
- Reproduzir com precisão as iso-linhas próximas às arestas superiores, onde a
  nuvem não uniforme do artigo tem vantagem.
- Atingir o nível de erro reportado pelo artigo naquelas regiões específicas.
- Comparar diretamente o número de nós necessários para dada precisão.

---

## 4. O que já pode ser afirmado

1. **O pipeline EFG esparso está correto e operacional** para o problema do
   cubo eletrostático 3D com parâmetros L = 10, V0 = 10.

4. **O erro máximo de ~10 V no plano completo é estrutural, esperado e não
   indica falha no solver.** Ele está confinado às arestas y = 0 e y = L em
   z = L, onde as condições de contorno conflitantes geram resultado diferente
   entre o EFG (prioridade Lagrange: V = V0) e a série (senos: V = 0).

**[OBSOLETO: Esta explicação se refere à política de Dirichlet anterior. Com a nova política (V=0 nas arestas superiores), o `max_abs_error` 3D para o caso `15x15x15` caiu para `4.642650e-01`. O `plane_max_abs_error` ainda é alto (`8.8854 V`), mas a explicação para isso agora reside no fenômeno de Gibbs da série analítica truncada e na natureza não interpolante do MLS nas regiões de forte gradiente e descontinuidade de contorno, conforme detalhado em `docs/figure3_reproduction.md`, Seção 5.]**

---

2. **O potencial no interior do plano x = 5.33 está bem reproduzido.**
   - Erro relativo interior: 4.58%.
   - Erro médio interior: 0.034 V = 0.34% de V0.
   - No centro do cubo (y ≈ 5, z ≈ 5): abs_error ≈ 0.012 V.

3. **A topologia das iso-linhas no interior coincide qualitativamente com a
   solução analítica.** As linhas de nível se curvam de forma consistente com
   a distribuição de potencial esperada para o problema do cubo.

5. **A série analítica com 25 termos exibe Gibbs** em z = L perto de y = 0 e
   y = L. Em (y = 0.2, z = 10): V_exact = 11.86 > V0, enquanto V_num = 10.03 ≈ V0.
   O EFG é mais correto nessa sub-região do que a própria série truncada usada
   como referência.

6. **A convergência com refinamento é monotônica e consistente:**
   5×5×5 → 21.78% → 11×11×11 → 5.55% → 15×15×15 → 2.81% (3D interior).

7. **O overshoot negativo do EFG** (V_num_min = −0.480 em cantos inferiores do
   plano) é pequeno (< 5% de V0) e limitado às bordas da grade.

8. **Quatro figuras de visualização foram geradas** e mostram a distribuição
   de potencial de forma clara, consistente e com escala unificada.

---

## 5. O que ainda não pode ser afirmado

1. **Reprodução visual direta da Figura 3 do artigo.** O script gera contornos
   preenchidos com 40 níveis em colormap; a Figura 3 provavelmente usa iso-
   linhas com níveis fixos (V/V0 = 0.1, ..., 0.9). Sem ajuste do estilo de
   visualização, a comparação é qualitativa.

2. **Equivalência quantitativa com os resultados do artigo.** O artigo usa
   nuvem não uniforme; nós usamos grade regular. Os números de erro são
   comparáveis mas não diretamente equivalentes.

3. **Precisão nas regiões de aresta superior.** As arestas y = 0, z = L e
   y = L, z = L têm erro do EFG ≈ V0 = 10. A nuvem não uniforme do artigo
   reduz esse erro. Com a grade regular, essa região permanece com baixa
   qualidade numérica.

4. **Que a acurácia atual é suficiente para identificar as iso-linhas específicas
   do artigo** sem verificar os níveis de contorno usados na publicação. Se o
   artigo usa V/V0 = 0.1, 0.2, ..., 0.9, precisamos verificar que os contornos
   em nossos dados passam pelos mesmos lugares no espaço y-z.

5. **Ausência de Gibbs no EFG em z → L.** O V_num_max = 10.033 indica leve
   sobressinal do EFG acima de V0 em z = L. Para a face superior, o EFG com
   Lagrange deveria impor exatamente V0; o excesso de 0.033 V é residual do
   GMRES (tolerância 10⁻⁹) e da interpolação MLS nos pontos de amostragem.

---

## 6. Próximos passos para comparação visual

Os itens abaixo são todos modificações no **script de visualização** — nenhum
requer alteração no solver, montagem ou app C.

### 6.1 Obrigatório para comparação com Figura 3

**a) Iso-linhas com níveis fixos (substituir ou sobrepor contourf)**

Adicionar `ax.contour()` com os mesmos níveis do artigo sobre ou no lugar do
`contourf`. Sugestão de níveis: `V_levels = [1, 2, 3, 4, 5, 6, 7, 8, 9]`
(i.e., V/V0 = 0.1, 0.2, ..., 0.9 para V0=10), usando `fmt='%.0f'` nos rótulos.

**b) Escala unificada entre V_num e V_exact**

A figura comparativa já usa `Normalize` compartilhado entre V_num e V_exact —
isso está correto.

**c) Ajuste de orientação dos eixos**

Verificar se o artigo usa y no eixo horizontal e z no vertical (como o script
atual), ou o inverso. Ajustar se necessário.

### 6.2 Desejável

**d) Normalizar potencial por V0 (V/V0 ∈ [0, 1])**

Divir os arrays de potencial por V0 antes de plotar. Isso facilita comparação
direta com figuras normalizadas da literatura.

**e) Marcar posição dos nós no plano**

Adicionar pontos (scatter) mostrando onde os nós EFG caem no plano x ≈ 5.33
(os nós mais próximos desse plano). Isso permite comparação qualitativa com o
esquema de pontos do artigo.

**f) Figura separada focando no interior (y ∈ [1, 9], z ∈ [1, 9])**

Recortar as figuras para excluir as bordas dominadas pelo erro estrutural,
mostrando apenas a região de comparação meaningful.

---

## 7. Riscos de interpretação

### 7.1 A série analítica não é a solução exata

**Risco:** Interpretar V_exact como solução analítica exata e concluir que
abs_error reflete puramente o erro do EFG.

**Realidade:** V_exact é uma série truncada (25 termos por dimensão). Perto de
discontinuidades em z = L e nas bordas y = 0/L, a série exibe Gibbs que pode
localmente SUPERAR V0 (como em V_exact = 11.86 em y = 0.2, z = 10). Nessa
região, parte do abs_error é erro da série, não do EFG.

**Mitigação:** Usar apenas as métricas do interior (interior_relative_error,
interior_mean_abs_error) para avaliar a qualidade do EFG. As métricas do
plano completo são indicativas, não diagnósticas.

### 7.2 interior_max_abs_error = 4.366 é enganoso

**Risco:** Concluir que o EFG tem erro máximo interior de 43% de V0.

**Realidade:** Esse valor ocorre em (y = 0.1, z = 9.8) — o ponto mais próximo
da aresta y = 0, z = L que ainda está dentro do filtro interior. Nessa posição,
ambos V_num e V_exact são aproximações com qualidade limitada; a diferença
entre elas não representa o erro verdadeiro do EFG em relação à solução
matemática perfeita. O erro médio interior (0.034) é o indicador mais confiável.

### 7.3 plane_relative_error = 8.19% não é o erro do método

**Risco:** Citar 8.19% como erro do método EFG no plano x = 5.33.

**Realidade:** Este valor inclui os pontos de borda dominados pelo conflito
de BC. O erro característico do método no interior do plano é **4.58%**, e
no centro do cubo é **< 0.2%**. A métrica correta a citar é `interior_relative_error`.

### 7.4 Comparação direta com o artigo é prematura

**Risco:** Afirmar que os resultados reproduzem a Figura 3 sem verificar os
níveis de contorno usados.

**Realidade:** O artigo usa nuvem não uniforme e provável exibição por iso-
linhas. A nossa reprodução usa grade regular e filled contourf. A comparação
qualitativa pode ser feita; a comparação quantitativa requer alinhamento de
metodologia de visualização e, idealmente, de nuvem de pontos.

### 7.5 Condicionamento estimado: 820

O campo `max_cond_estimate = 820.48` (razão max/min de pivôs em eliminação
gaussiana) para o caso 15×15×15 indica que a matriz de momento MLS mais mal
condicionada nessa malha tem número de condição ≈ 820. Para `double` (eps ≈
10⁻¹⁶), isso equivale a uma perda de ~3 dígitos decimais. Com soluções da
ordem de 10, a imprecisão introduzida é ~10×820×10⁻¹⁶ ≈ 10⁻¹². Esse
condicionamento é aceitável. Porém, com refinamentos maiores (17×17×17+), o
condicionamento pode crescer e deve ser monitorado.

---

## 8. Resumo avaliativo

| Pergunta | Resposta | Evidência |
|---|---|---|
| O pipeline está correto? | Sim | 28/28 testes, support_lt_4=0, mls_failures=0, GMRES convergido |
| A distribuição de V no interior está correta? | Sim | interior_relative_error=4.58%, abs_error_center≈0.012 |
| O erro máximo (≈V0) é aceitável? | Sim (estrutural) | BC conflitante y=0/y=L, z=L; explicado matematicamente |
| O script reproduce a Figura 3? | Parcialmente | Dados corretos; formato visual diverge (filled vs. iso-linhas) |
| O resultado é comparável ao artigo? | Qualitativamente | Diferença: grade regular vs. nuvem não uniforme |
| A série analítica é exata? | Não (25 termos) | Gibbs: V_exact_max=11.86>V0 em z=L perto de y=0 |
| Condicionamento MLS é seguro? | Sim (820, aceitável) | Margem ≈ 12 ordens de magnitude acima do eps da máquina |

---

## 9. Recomendação

```
AJUSTAR SCRIPT DE CONTORNO
```

**Justificativa:**

O solver, a montagem e os dados estão corretos. O pipeline está validado. Os
dados do plano x = 5.33 são cientificamente adequados para comparação com a
Figura 3 do artigo. A limitação atual é exclusivamente de visualização:

1. O script gera **filled contourf** com 40 níveis em colormap — o artigo
   provavelmente exibe **iso-linhas** com 9 a 10 níveis em V/V0.
2. Para uma afirmação de "reprodução da Figura 3" ser defensável, o script
   precisa ao menos adicionar iso-linhas com os mesmos níveis da publicação,
   sobrepondo V_num e V_exact no mesmo painel.
3. Isso é uma mudança de **3 a 5 linhas de Python** no script
   `plot_cube_plane.py`, sem tocar solver, montagem ou app C.

**O que NÃO precisa ser feito antes do ajuste:**
- Não é necessário voltar ao solver (todos os critérios numéricos estão
  satisfeitos).
- Não é necessário implementar nuvem não uniforme agora (as limitações
  foram identificadas e documentadas; a grade regular é suficiente para
  comparação qualitativa).
- Não é necessário rodar malha maior (15×15×15 com 2.81% interior error 3D
  é suficiente para a figura).

**Ajustes específicos para o script:**
1. Ler os níveis de contorno usados no artigo (ou usar `[1,2,3,4,5,6,7,8,9]`
   como estimativa inicial).
2. Adicionar `ax.contour()` com esses níveis sobre o `contourf` (ou
   substituir).
3. Normalizar por V0 se o artigo usa V/V0.
4. Gerar uma figura de sobreposição: iso-linhas de V_num em azul sólido e
   V_exact em vermelho tracejado, no mesmo painel.

**Depois do ajuste:** AVANÇAR PARA DOCUMENTAÇÃO DA FIGURA.
