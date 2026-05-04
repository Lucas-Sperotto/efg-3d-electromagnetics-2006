# Parâmetros numéricos de referência

Este arquivo consolida os valores numéricos citados na tradução do artigo para orientar a implementação didática em C.

Os valores abaixo são uma **referência inicial de reprodução**. Eles não devem ser tratados como validação automática por si só. Uma validação numérica futura ainda deve comparar a solução implementada com a solução analítica, com os dados do artigo ou com os gráficos correspondentes.

---

## Parâmetros do problema de referência

| Parâmetro | Valor no artigo | Fonte na tradução | Constante C sugerida |
|---|---:|---|---|
| Lado do cubo | `10` | `docs/03_resultados_numericos.md`, seção 3.1 | `EFG_CUBE_L` |
| Potencial imposto na parede superior | `10 V` | `docs/03_resultados_numericos.md`, seção 3.1 | `EFG_CUBE_V0` |
| Plano de corte para a Fig. 3 | `x = 5.33` | `docs/03_resultados_numericos.md`, seção 3.3 | `EFG_SLICE_X` |
| Células de integração | `15 x 15 x 15` | `docs/03_resultados_numericos.md`, seção 3.3 | `EFG_GAUSS_CELLS_X`, `EFG_GAUSS_CELLS_Y`, `EFG_GAUSS_CELLS_Z` |
| Quadratura de Gauss | ordem `2` | `docs/03_resultados_numericos.md`, seção 3.3 | `EFG_GAUSS_ORDER` |
| Função peso principal | spline cúbica | `docs/02_formulacao_matematica.md`, seção 2.2; `docs/03_resultados_numericos.md`, seção 3.3 | `EFG_WEIGHT_CUBIC_SPLINE` |
| Função peso de comparação | quadrática | `docs/02_formulacao_matematica.md`, seção 2.2; `docs/03_resultados_numericos.md`, seção 3.4 | `EFG_WEIGHT_QUADRATIC` |
| Domínio de influência usado nos resultados | cúbico por produto tensorial | `docs/02_formulacao_matematica.md`, seção 2.2; `docs/03_resultados_numericos.md`, seção 3.3 | `EFG_INFLUENCE_CUBIC_TENSOR` |
| Número de vizinhos mais próximos | `k = 3` | `docs/03_resultados_numericos.md`, seção 3.2 | `EFG_NEIGHBOR_K` |
| Fator de suporte | `1.5` vezes a distância ao terceiro vizinho | `docs/03_resultados_numericos.md`, seção 3.2 | `EFG_SUPPORT_SCALE` |
| Erro relativo EFG de referência | `1.40%` | `docs/03_resultados_numericos.md`, seção 3.3; `docs/04_conclusao_referencias.md` | `EFG_REFERENCE_ERROR_EFG` |
| Erro relativo FEM de referência | `3.56%` | `docs/03_resultados_numericos.md`, seção 3.3; `docs/04_conclusao_referencias.md` | `EFG_REFERENCE_ERROR_FEM` |

---

## Convenção sugerida para valores escalares

| Constante C sugerida | Valor inicial sugerido | Observação |
|---|---:|---|
| `EFG_CUBE_L` | `10.0` | Lado do cubo do problema analítico. |
| `EFG_CUBE_V0` | `10.0` | Potencial prescrito na face `z = L`, em volts. |
| `EFG_SLICE_X` | `5.33` | Plano usado para reproduzir a fatia `yz` da Fig. 3. |
| `EFG_GAUSS_CELLS_X` | `15` | Número de células de integração na direção `x`. |
| `EFG_GAUSS_CELLS_Y` | `15` | Número de células de integração na direção `y`. |
| `EFG_GAUSS_CELLS_Z` | `15` | Número de células de integração na direção `z`. |
| `EFG_GAUSS_ORDER` | `2` | Ordem da quadratura por direção. |
| `EFG_NEIGHBOR_K` | `3` | Terceiro vizinho mais próximo. |
| `EFG_SUPPORT_SCALE` | `1.5` | Fator multiplicativo do domínio de influência. |
| `EFG_REFERENCE_ERROR_EFG` | `0.0140` | Valor em escala fracionária, equivalente a `1.40%`. |
| `EFG_REFERENCE_ERROR_FEM` | `0.0356` | Valor em escala fracionária, equivalente a `3.56%`. |

---

## Observações para implementação

- A função peso spline cúbica é a referência principal dos resultados numéricos do artigo.
- A função peso quadrática aparece como comparação e não substitui automaticamente a escolha principal.
- O domínio de influência cúbico usa produto tensorial das funções peso em cada direção.
- O fator `1.5` é aplicado à distância até o terceiro vizinho mais próximo.
- TODO: confirmar manualmente se a distância até o terceiro vizinho deve ser sempre euclidiana no artigo original.
- TODO: documentar separadamente qualquer alteração desses valores usada em experimentos próprios.
