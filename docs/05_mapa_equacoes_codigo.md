# Mapa entre equações e código C

Este mapa orienta a implementação didática do artigo. Ele relaciona as equações e conceitos da tradução com arquivos e funções C sugeridos.

TODO: a numeração das equações ainda precisa ser revisada manualmente contra o artigo original. Os arquivos Markdown atuais citam algumas equações pelo número do artigo, mas nem todas as equações estão rotuladas na tradução.

---

## Tabela de mapeamento

| Equação do artigo | Significado | Arquivo C sugerido | Função C sugerida | Status |
|---|---|---|---|---|
| Eq. (1) | Forma forte de Poisson para o potencial eletrostático | `src/assembly.c` | `assembly_poisson_strong_form_reference()` | Documentação pendente; não implementar diretamente antes da montagem. |
| Eqs. (2)-(3) | Condições de contorno de Dirichlet e Neumann homogênea | `src/assembly.c` | `assembly_apply_boundary_model()` | TODO; depende da modelagem da fronteira. |
| Eq. (4) | Funcional com multiplicadores de Lagrange | `src/assembly.c` | `assembly_lagrange_functional_reference()` | TODO; funcional é referência conceitual, não necessariamente função executável. |
| Eqs. (5)-(6) | Forma fraca com funções teste para `V` e `lambda` | `src/assembly.c` | `assembly_weak_form()` | TODO; depende de quadratura, MLS e dados de fronteira. |
| Eq. (7) | Sistema matricial com blocos `K`, `G`, `G^T` e multiplicadores | `src/assembly.c` | `assembly_build_lagrange_system()` | TODO; depende de armazenamento de matrizes e vetores. |
| Eq. (8) | Matriz de rigidez `K` | `src/assembly.c` | `assembly_stiffness_matrix()` | TODO; depende das derivadas das funções de forma. |
| Eq. (9) | Blocos de fronteira `G`, vetor `b` e vetor de fonte `f` | `src/assembly.c` | `assembly_boundary_terms()` | TODO; depende da fronteira de Dirichlet e da regra para `N_i`. |
| Eq. (10) | Função de forma MLS `Phi_I(x)` | `src/mls.c` | `mls_shape_function()` | TODO; não implementar antes de definir suporte e álgebra local. |
| Eq. (10), conforme citado em `docs/03_resultados_numericos.md` | Função peso spline cúbica | `src/weight.c` | `weight_cubic_spline()` | Implementado. TODO: revisar numeração porque há conflito com a numeração provável de MLS. |
| Eq. (11), conforme citado em `docs/03_resultados_numericos.md` | Função peso quadrática | `src/weight.c` | `weight_quadratic()` | Implementado. TODO: revisar numeração contra o artigo original. |
| Eqs. (13)-(15) | Solução analítica do cubo por série dupla | `src/analytical.c` | `analytical_potential_cube()` | Implementado com truncamento finito. |
| Eq. (16) | Norma relativa do erro por integral de domínio | `src/error_norm.c` | `relative_error_norm_domain_integral()` / `relative_error_norm_domain_integral_with_order()` | Implementado com Gauss 2x2x2 como padrão e ordem 1..8 explícita para estudos de sensibilidade. |
| Referência [7] | Solução do sistema não definido positivo por GMRES | `src/gmres.c` | `gmres_solve()` | TODO; não implementar antes da interface de matriz/vetor. |

---

## Ordem recomendada para implementação futura

1. Revisar a numeração de equações contra o artigo original.
2. Definir estruturas de nós, fronteiras e valores prescritos.
3. Definir suporte nodal e domínio de influência cúbico.
4. Implementar MLS e derivadas necessárias.
5. Implementar quadratura e montagem de `K`, `G`, `b` e `f`.
6. Implementar ou adaptar GMRES didático.
7. Implementar norma relativa de erro e comparar com a solução analítica.

---

## Notas de consistência

- Toda função C que implementar uma equação deve citar o número da equação, o arquivo Markdown de referência e o significado físico ou matemático.
- Onde a numeração ainda estiver ambígua, o comentário deve incluir `TODO: confirmar número da equação no artigo original`.
- Não usar este mapa como substituto do artigo: ele é uma camada de organização para reduzir erro de implementação.
