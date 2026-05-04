# Notas sobre a ESFLib — referência técnica para implementação do MLS/EFG em C

---

## 1. O que é a ESFLib

A **ESFLib** (*EFG Shape-Function Library*) é uma biblioteca de referência para o cálculo das funções de forma do método Element-Free Galerkin (EFG), escrita em **ANSI C** segundo princípios de orientação a objetos.

- **Autor:** Petr Krysl, Computational Mechanics Group, Northwestern University
- **Ano:** 1997
- **Copyright:** © Petr Krysl, 1997
- **Localização no repositório:** `docs/ref/ESFLib/` (extraída de `docs/ref/ESFLib.tar.gz`)
- **Referência original:** P. Krysl e T. Belytschko, "ESFLIB: A library to compute the element free Galerkin shape functions", submetido a *Computer Methods in Applied Mechanics and Engineering*, 1997.

A biblioteca é **independente de estrutura de dados**: toda conectividade nodal é fornecida via callbacks, o que a torna altamente configurável sem exigir estruturas específicas do aplicativo chamador.

---

## 2. Arquivos relevantes para consulta

| Arquivo | Conteúdo |
|---|---|
| `doc/esflib.html` | Documentação principal: design, API, exemplos de uso |
| `src/sf.c` | Implementação central do MLS: montagem de A(x), solução do sistema, cálculo de Φ_I e suas derivadas |
| `src/wf.c` | Todas as funções peso disponíveis (quartic spline, quadrática, Gaussiana, produto tensorial) com derivadas de 1ª e 2ª ordem |
| `src/pb3.c` | Base polinomial 3D: vetores p(x) e dp/dx para ordens constante, linear, quadrática, cúbica e quártica |
| `src/esflibpub.h` | API pública: tipos, enumerações de erro, protótipos de todas as funções |
| `src/esflibtest.c` | Programa de teste: partição da unidade, soma das derivadas, testes de desempenho |

---

## 3. O projeto atual NÃO copiará o código da ESFLib

A implementação didática deste repositório será escrita de forma **independente**, em C11, específica para o problema do artigo de Parreira et al. (2006). A ESFLib **não será compilada nem vinculada** ao código deste projeto.

Razões:
1. O objetivo é didático: cada função deve ser rastreável a uma equação numerada da `docs/`.
2. A ESFLib é genérica (1D/2D/3D) e usa macros de dimensão em tempo de compilação, o que dificulta a leitura para fins de aprendizado.
3. A licença não está explicitamente declarada como permissiva; redistribuição requer verificação (ver §5).

---

## 4. O que a ESFLib orienta na implementação

### 4.1. Organização do cálculo das funções de forma

`sf.c` revela a sequência canônica para calcular Φ_I(x) (cf. equações (10)–(14) do artigo):

1. Buscar os nós I cujo suporte contém x (callback `find_overlaps`).
2. Para cada nó I encontrado: calcular w(r_I) e p(x_I).
3. Montar A(x) = Σ w_I p(x_I) p^T(x_I) — matriz 4×4 para base linear 3D.
4. Resolver A(x) a(x) = B_I(x) por LU com pivotamento parcial (padrão) ou QR.
5. Calcular Φ_I(x) = p^T(x) A^{-1}(x) B_I(x).
6. Para derivadas: diferenciar a relação acima, propagando dA/dx e dB_I/dx.

A biblioteca oferece três decomposições para A: `ESFLIB_LU_DECOMP` (padrão, ~14% mais rápida), `ESFLIB_QR_DECOMP` (mais robusta) e `ESFLIB_CH_DECOMP` (Cholesky, para A simétrica definida positiva).

### 4.2. Testes de partição da unidade

`esflibtest.c` verifica que:

```
Σ_I Φ_I(x) = 1   para todo x
```

Este teste deve ser replicado na implementação própria antes de qualquer integração numérica.

### 4.3. Testes de soma das derivadas

`esflibtest.c` também verifica:

```
Σ_I ∂Φ_I/∂x = 0
Σ_I ∂Φ_I/∂y = 0
Σ_I ∂Φ_I/∂z = 0
```

Essas condições decorrem da consistência linear do espaço de aproximação e são necessárias para garantir que a matriz de rigidez K seja construída corretamente.

### 4.4. Estrutura de vizinhança (domínio de influência)

A ESFLib delega a busca de vizinhos ao aplicativo via três callbacks:

```c
find_overlaps(SS, x, y, z)      /* localiza todos os nós com suporte em x */
first_overlapping(SS)            /* retorna o primeiro nó encontrado */
next_overlapping(SS)             /* itera pelos demais */
```

Para a implementação própria: o critério de inclusão para domínio cúbico é `|x−x_I| ≤ d_mI` **em cada componente** separadamente (não pela norma euclidiana). A ESFLib confirma isso ao passar `rx`, `ry`, `rz` independentemente para as funções peso.

### 4.5. Cuidados com condicionamento da matriz A

A ESFLib expõe:
- `ESFLIB_A_cond_est()` — estimativa do número de condição de A(x).
- Código de erro `ESFLIB_FAILED_SOLVE` — A singular ou mal-condicionada.
- Código de erro `ESFLIB_BAD_PREC` — solução obtida com possível perda de precisão.
- Parâmetro `ESFLIB_set_min_num_neighbors()` — número mínimo de nós no suporte.

Na implementação própria: verificar que A seja bem condicionada é **obrigatório** antes de chamar o solver. Para base linear 3D com 4 termos, o suporte deve conter **pelo menos 4 nós** para que A seja não-singular (na prática, usar ≥ 6–8 nós é mais seguro).

---

## 5. Alerta: licença e redistribuição

> **Verificar a licença antes de qualquer redistribuição.**

A ESFLib exibe `(C) Petr Krysl, 1997` na documentação (`doc/esflib.html`). Nenhuma licença explícita (MIT, GPL, LGPL, BSD) foi encontrada no arquivo compactado. Isso significa que, pelo padrão de direitos autorais, **todos os direitos são reservados** a menos que o autor conceda permissão expressa.

Consequências práticas para este projeto:
- O arquivo `docs/ref/ESFLib.tar.gz` e a pasta `docs/ref/ESFLib/` devem ser tratados como **material de referência de uso pessoal**, não como código de terceiros redistribuível.
- **Não incorporar** trechos de código da ESFLib diretamente nos arquivos `.c` ou `.h` deste repositório.
- Se o repositório for tornado público, avaliar se `docs/ref/ESFLib/` deve ser excluído do versionamento (via `.gitignore`) ou mantido apenas localmente.

---

## 6. Correspondência entre funções da ESFLib e o artigo de 2006

A tabela abaixo identifica quais funções da ESFLib correspondem às equações do artigo (quando aplicável):

| Símbolo no artigo | Descrição | Equivalente na ESFLib | Observação |
|---|---|---|---|
| w(r) spline cúbica | Eq. (10) do artigo | **sem equivalente direto** | ESFLib tem `ESFLIB_qs_w` (quartic spline: 1−6r²+8r³−3r⁴), não a cúbica do artigo |
| w(r) quadrática | Eq. (11) do artigo | `ESFLIB_sqtp_w` (3D tensor product) | Fórmula idêntica: 1−2r² e 2(1−r)² |
| w(r_x)w(r_y)w(r_z) | Produto tensorial | Implementação interna de todos `_tp_` | Confirmado em `wf.c` |
| p^T(x) = [1,x,y,z] | Base linear 3D | `ESFLIB_poly_basis(..., 4, ...)` | Ver `pb3.c`, nterms=4 |
| A(x) | Matriz momento | Montada em `sf.c` | LU com pivotamento por padrão |
| Φ_I(x) = p^T A^{-1} B_I | Função de forma | `ESFLIB_shape_func_at()` | API pública em `esflibpub.h` |
| ∂Φ_I/∂x | Derivadas | `ESFLIB_inspect_row(..., ESFLIB_ORDER_1, ESFLIB_X, ...)` | Mesma chamada para y, z |

**Atenção:** A função peso spline cúbica usada no artigo (Eq. 10) **não existe na ESFLib**. A ESFLib tem a quartic spline `1 − 6r² + 8r³ − 3r⁴` (diferente) e a quadrática composta (que corresponde à Eq. 11). A spline cúbica do artigo deve ser implementada do zero na função `weight_cubic_spline()` do projeto.

---

**Anterior:** [6. Condições de contorno e multiplicadores de Lagrange](06_condicoes_de_contorno_e_lagrange.md) | **Voltar ao:** [README](../README.md)
