# Revisão da `docs/` como baseline científica para implementação em C

---

## 1. Numeração de equações

**INCONSISTÊNCIA CRÍTICA — as equações nos `.md` não têm numeração.**

O texto de `03_resultados_numericos.md` referencia explicitamente "Eq. (7)", "Eqs. (8) e (9)", "Eq. (10)", "Eq. (11)" e "Eq. (16)", mas nenhuma equação em nenhum dos arquivos `.md` possui etiqueta `(n)`.

O `AGENTS.md` (regra 3) exige que *toda* função C cite o número da equação. Essa regra é **impossível de seguir** enquanto as equações não estiverem numeradas na documentação.

**Recomendação:** Adicionar rótulos `(1)`, `(2)`, … a todas as equações de `02_formulacao_matematica.md` e `03_resultados_numericos.md`, preservando a numeração do artigo original (deduzível pelas referências internas). Mapeamento provável:

| Eq. no texto | Equação |
|---|---|
| (1) | Poisson: ∇·(ε∇V) = −ρ |
| (2)–(3) | Condições de contorno |
| (4) | Funcional Π(V, λ) |
| (5)–(6) | Forma fraca (variacional) |
| (7) | Sistema matricial [K G^T; G 0] |
| (8) | K_ij |
| (9) | G_ij, b_i, f_i |
| (10) | Φ_I(x) = p^T A^{-1} B_I |
| (10)–(11)? | w(r) spline cúbica / quadrática |
| (16) | Norma relativa do erro |

A numeração exata precisa ser reconstituída confrontando com o PDF original (DOI: 10.1109/TMAG.2006.872014).

---

## 2. Ambiguidades que podem prejudicar a implementação

**A. Simplificação pela função delta de Dirac — sem fórmula resultante** (`03_resultados_numericos.md`, §3.3)

> "duas integrais de superfície da Eq. (7) não precisam ser realizadas"

Não está escrito qual é a expressão resultante. O implementador precisa deduzir:
- G_ij = Φ_j(x_i) (avaliação pontual da função de forma no nó i de fronteira)
- b_i = V_d(x_i) (valor prescrito no nó i)

**Risco:** Implementação errada da imposição de contorno de Dirichlet — o erro mais comum em EFG.
**Recomendação:** Adicionar subseção com as formas simplificadas explícitas.

**B. Parâmetro d_mI: isotrópico ou anisotrópico?** (`02_formulacao_matematica.md`)

As fórmulas `r_x = |x−x_I|/d_mI`, `r_y = |y−y_I|/d_mI`, `r_z = |z−z_I|/d_mI` usam um único escalar `d_mI` para as três direções. No texto (`03_resultados_numericos.md`) define-se `d_mI = 1,5 × distância ao 3º vizinho`. Mas que distância? Euclidiana, máxima das componentes? Para nós na fronteira (onde a distribuição é 2D), a distância 3D pode subestimar o suporte.

**Recomendação:** Declarar explicitamente que `d_mI` é a distância euclidiana ao k-ésimo vizinho.

**C. Domínio de influência esférico vs. cúbico — definição de r não fornecida para o caso esférico** (`02_formulacao_matematica.md`)

O texto apresenta ambas as opções mas não define `r` para o caso esférico. Para o caso cúbico usa-se o produto tensorial `w(r_x)w(r_y)w(r_z)`. Para o esférico, presume-se `r = ‖x − x_I‖/d_mI`, mas isso não está escrito.

**Recomendação:** Adicionar a definição `r = √(r_x² + r_y² + r_z²)` para o domínio esférico.

**D. n_x em A(x): critério de inclusão não formalizado** (`02_formulacao_matematica.md`)

```
A(x) = Σ_{I=1}^{n_x} w(x−x_I) p(x_I) p^T(x_I)
```

`n_x` é o número de nós cujo suporte cobre `x`. O critério de inclusão (`‖x − x_I‖ ≤ d_mI`? ou por componente?) não está escrito. Para o domínio cúbico o critério muda. Sem a definição explícita, o loop de busca nodal é ambíguo.

**E. Que equação é referenciada como "Eq. (7)" em §3.3?**

O texto diz "duas integrais de superfície da Eq. (7)". Ao remontar a numeração, Eq. (7) parece ser o sistema matricial. Mas o sistema matricial contém submatrizes, não "integrais de superfície" isoladas. A referência provavelmente é às definições de G_ij e b_i. Precisa ser esclarecido.

---

## 3. Símbolos sem definição explícita

| Símbolo | Onde aparece | Problema |
|---|---|---|
| **x** (duplo uso) | `Φ_I(x) = …` e `p^T(x) = [1, x, y, z]` | `x` é vetor posição e também coordenada escalar — notação sobrecarregada |
| **N_λ** | `λ(x) ≈ Σ_{i=1}^{N_λ} N_i λ_i` | Nunca definido. Implicitamente = nº de nós de fronteira Dirichlet quando se usa delta de Dirac |
| **N_i(x)** | Interpolante de λ | Natureza não definida até §3.3 (descobrimos que é delta de Dirac). No §2 parece uma função arbitrária |
| **B_I(x)** | `B_I(x) = w(x−x_I) p(x_I)` | Dimensão não declarada. É um vetor coluna 4×1 (mesmo tamanho de p) |
| **A^{-1}(x)** | Eq. Φ_I | Condicionamento/singularidade mencionado verbalmente mas sem critério numérico |
| **w(x−x_I)** | Soma em A(x) | Argumento vetorial para função definida com argumento escalar `r`. Não é dito que `r = ‖x−x_I‖/d_mI` |
| **ρ** no benchmark | `f_i = ∫ ρ Φ_i dΩ` | No problema de referência é equação de Laplace (ρ = 0), mas isso não é declarado em §3 |

---

## 4. Solução analítica do cubo — completude

**A. Condições de contorno não escritas formalmente**

A solução é dada diretamente, mas as BCs não são listadas:
- V = 0 em x=0, x=L, y=0, y=L, z=0
- V = V_0 em z=L

A descrição verbal diz "parede superior com potencial diferente" mas não formaliza quais paredes são aterradas.

**Risco:** Implementação com orientação de eixo errada (confundir "superior" com z vs y).

**B. Truncamento da série infinita — sem critério**

A série é dupla e infinita (n, m ímpares). Nenhum critério de truncamento é fornecido. Para fins de implementação do verificador analítico em C, é necessário saber quantos termos usar.

**Recomendação:** Adicionar N_max = 20–50 (n, m = 1, 3, 5, …, N_max) com nota sobre a convergência rápida para pontos interiores e lenta para cantos (onde há singularidade de aresta).

**C. Origem do coeficiente 160 não explicada**

`a_mn = 160/(mn π² sinh(10 k_mn))` — o 160 é `16 × V_0 = 16 × 10`. A derivação (projeção na base de Fourier dupla) não está documentada. Para V_0 e L genéricos a fórmula seria `16 V_0/(mn π² sinh(k_mn L))`, o que é importante para testar com outros parâmetros.

**D. Plano x = 5,33 da Fig. 3 — sem justificativa**

Para um cubo de lado L=10 o plano central seria x=5. O valor 5,33 não é explicado.

**Risco baixo**, mas impede reproduzir a Fig. 3 exata sem saber a justificativa.

---

## 5. Parâmetros numéricos dispersos que devem ir para constantes no código

Estão espalhados em dois arquivos e texto corrido, sem consolidação:

| Parâmetro | Valor | Seção fonte |
|---|---|---|
| Dimensão do cubo L | 10 | §3.1 |
| Potencial prescrito V_0 | 10 V | §3.1 |
| Número de células de Gauss | 15×15×15 | §3.3 |
| Ordem de Gauss | 2 (por dimensão) | §3.3 |
| k-vizinhos mais próximos | 3 | §3.2 |
| Fator d_mI | 1,5 × dist. ao 3º vizinho | §3.2 |
| Tipo de função peso principal | spline cúbica | §3.3 |
| Tipo de domínio de influência | cúbico (produto tensorial) | §3.3 |
| Plano de corte para Fig. 3 | x = 5,33 | §3.3 |
| Erro EFG esperado (validação) | 1,40% | §3.3 / §4 |
| Erro FEM (referência) | 3,56% | §3.3 / §4 |

**Recomendação:** Criar um arquivo `docs/parametros_numericos.md` (ou tabela em `03_resultados_numericos.md`) consolidando todos esses valores, com os nomes de constantes sugeridos para o código C (ex: `#define EFG_L 10.0`, `#define EFG_V0 10.0`, etc.).

---

## 6. Arquivo de mapeamento equação → código

**Não existe.** O `AGENTS.md` (regra 3) exige que toda função C cite número de equação e arquivo `.md`, mas:
1. As equações não têm números (problema 1 acima).
2. Não existe nenhum arquivo `docs/mapa_equacoes.md` ou similar.

**Recomendação:** Criar `docs/mapa_equacoes.md` com tabela mínima:

| Eq. | Expressão | Função C sugerida | Arquivo `.md` |
|---|---|---|---|
| (4) | Π(V,λ) | — (funcional, não implementado diretamente) | `02_formulacao_matematica.md` |
| (7) | Sistema [K G^T; G 0] | `build_system()` | `02_formulacao_matematica.md` |
| (8) | K_ij | `stiffness_matrix()` | `02_formulacao_matematica.md` |
| (10) | Φ_I = p^T A^{-1} B_I | `efg_shape_function()` | `02_formulacao_matematica.md` |
| (10) | w(r) spline cúbica | `weight_cubic_spline()` | `02_formulacao_matematica.md` |
| (16) | norma relativa do erro | `relative_error_norm()` | `03_resultados_numericos.md` |
| — | Solução analítica | `analytic_solution()` | `03_resultados_numericos.md` |

---

## Resumo de prioridades

| # | Problema | Impacto | Ação |
|---|---|---|---|
| P1 | Equações sem numeração | **Bloqueador** — inviabiliza regra 3 do AGENTS.md | Numerar todas as equações |
| P2 | Simplificação delta de Dirac sem expressão resultante | **Alto** — erro direto na imposição de BC | Documentar G_ij e b_i simplificados |
| P3 | BCs da solução analítica não formalizadas | **Alto** — orientação de eixo errada | Listar as 6 faces e seus potenciais |
| P4 | Arquivo de mapa equação → código inexistente | **Alto** — inviabiliza AGENTS.md | Criar `docs/mapa_equacoes.md` |
| P5 | N_λ, N_i, x (duplo uso) indefinidos | **Médio** — confusão de dimensões | Adicionar glossário de símbolos |
| P6 | Parâmetros numéricos dispersos | **Médio** — risco de valores inconsistentes no código | Consolidar em tabela única |
| P7 | Truncamento da série analítica sem critério | **Médio** — implementação do verificador incompleta | Adicionar N_max recomendado |
| P8 | d_mI: distância euclidiana não declarada | **Baixo-médio** — erro silencioso em fronteiras | Explicitar a norma usada |
| P9 | Plano x=5,33 sem justificativa | **Baixo** — impede reproduzir Fig. 3 exata | Documentar ou corrigir |
| P10 | Coeficiente 160 sem derivação | **Baixo** — implementação correta mas frágil | Adicionar fórmula geral em V_0, L |
