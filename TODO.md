# TODO — Fechamento científico e computacional do EFG 3D

Este arquivo organiza o estado atual e o plano de continuidade do projeto `efg-3d-electromagnetics-2006`.

O objetivo principal é reproduzir, de forma didática e verificável, o problema eletrostático tridimensional do artigo:

> The Element-Free Galerkin Method in Three-Dimensional Electromagnetic Problems  
> IEEE Transactions on Magnetics, 2006.

A prioridade agora não é criar investigações secundárias, mas consolidar a formulação numérica para permitir refinamento suficiente e reprodução confiável das figuras do artigo.

---

## 0. Estado atual do projeto

### 0.1. Base já implementada

- [x] Estrutura C11 com CMake.
- [x] Funções peso escalares.
- [x] Derivadas escalares das funções peso.
- [x] Peso tensorial 3D com gradiente físico.
- [x] Cálculo de suporte nodal baseado no terceiro vizinho multiplicado por 1,5.
- [x] Solver local 4x4 para MLS linear 3D.
- [x] MLS linear 3D sem derivadas.
- [x] MLS linear 3D com gradientes.
- [x] Quadratura de Gauss 2x2x2.
- [x] Montagem local da matriz de rigidez.
- [x] Montagem global densa da matriz de rigidez.
- [x] Montagem densa das restrições de Dirichlet via multiplicadores de Lagrange.
- [x] Montagem densa do sistema aumentado:

```text
  [ K  G^T ]
  [ G   0  ]
```

- [x] Solver denso didático por eliminação de Gauss.
- [x] Teste integrado artificial com solução linear.
- [x] Gerador de problema cúbico regular.
- [x] Geração de Dirichlet a partir de grade regular.
- [x] Geração de Dirichlet a partir de nuvem arbitrária de nós.
- [x] Aplicação `reproduce_cube_smoke`.
- [x] Script de análise de erro por região.
- [x] Script de estudo de refinamento regular.
- [x] Script de plotagem da convergência.
- [x] Primeira tentativa de nuvem não uniforme refinada no topo.

---

## 1. Diagnóstico importante

A primeira tentativa de nuvem não uniforme não melhorou o resultado.

Comparação observada:

```text
Regular 11x11x11 com integração 15x15x15:
nodes: 1331
constraints: 602
internal relative error: ~0.05294

Não uniforme base_n=7 refine_n=11:
nodes: 1308
constraints: 786
internal relative error: ~0.15300
```

Conclusão provisória:

- a nuvem não uniforme ficou dominada por pontos de contorno;
- houve aumento relativo das restrições de Dirichlet;
- a distribuição pode ter piorado o condicionamento local;
- ainda não há diagnóstico de cardinalidade e condicionamento do MLS;
- ainda não há estrutura esparsa para resolver casos suficientemente refinados.

Portanto, a linha de trabalho com nuvens não uniformes deve ficar pausada.

---

## 2. Correção de rumo

Antes de continuar refinando a geometria ou tentando reproduzir a figura final, é necessário garantir:

1. cada ponto de integração possui número suficiente de nós ativos no suporte;
2. a cardinalidade local não é inferior ao número de termos da base;
3. a matriz de momento MLS não fica singular ou mal condicionada;
4. a montagem global deve migrar de densa para esparsa;
5. o sistema aumentado deve ser resolvido por método iterativo adequado, inicialmente GMRES;
6. só depois disso deve-se buscar refinamento suficiente para reproduzir a figura do artigo.

Para a base linear 3D:

```text
p(x,y,z) = [1, x, y, z]
```

a cardinalidade mínima teórica é:

```text
4 nós ativos
```

mas, na prática, é recomendável ter mais do que 4 nós no suporte para reduzir singularidade e mau condicionamento.

---

## 3. Prioridade 1 — Diagnóstico de conectividade MLS

### 3.1. Objetivo

Criar ferramentas para avaliar a conectividade local de cada ponto de integração e de cada ponto de avaliação.

Para cada ponto analisado, calcular:

- número de nós ativos no suporte;
- se a cardinalidade é menor que 4;
- se a cardinalidade é pequena, por exemplo menor que 8;
- se a matriz de momento MLS é singular ou quase singular;
- distribuição estatística da conectividade.

### 3.2. Tarefas

- [x] Criar módulo de diagnóstico MLS.
- [x] Criar função para contar nós ativos no suporte de um ponto.
- [x] Criar função para avaliar estatísticas de conectividade em uma grade de pontos.
- [x] Criar teste unitário para conectividade em grade regular.
- [ ] Criar teste unitário para conectividade em nuvem não uniforme.
- [ ] Criar app ou script que gere CSV com:

  ```text
  x,y,z,active_nodes,status
  ```

- [ ] Gerar relatório para os casos:

  ```text
  3x3x3
  5x5x5
  7x7x7
  9x9x9
  11x11x11
  ```

### 3.3. Critério de parada

Não avançar para refinamentos maiores enquanto houver pontos de integração ou amostragem com:

```text
active_nodes < 4
```

Também registrar pontos com:

```text
active_nodes < 8
```

como região numericamente suspeita.

---

## 4. Prioridade 2 — Diagnóstico da matriz de momento MLS

### 4.1. Objetivo

Avaliar a qualidade da matriz:

```text
A(x) = Σ_I w_I p_I p_I^T
```

em pontos de integração e amostragem.

### 4.2. Tarefas

- [x] Criar função para montar `A(x)` sem calcular necessariamente as funções de forma.
- [x] Criar função para estimar condicionamento ou pelo menos detectar pivôs pequenos durante a eliminação.
- [ ] Registrar:

  - menor pivô;
  - falha de solução;
  - número de nós ativos;
  - coordenada do ponto;
  - tipo de ponto: integração, avaliação ou Dirichlet.
- [ ] Gerar CSV diagnóstico.
- [ ] Comparar grade regular e nuvem não uniforme.

### 4.3. Critério de parada

Não buscar reprodução quantitativa enquanto houver falhas frequentes em MLS ou regiões com matriz de momento quase singular.

---

## 5. Prioridade 3 — Matriz esparsa

### 5.1. Justificativa

A montagem densa permitiu validar o pipeline, mas não permite refinamento suficiente.

Exemplo:

```text
15x15x15 nós = 3375 nós
restrições de superfície ≈ 1178
sistema aumentado ≈ 4553 x 4553
```

Uma matriz densa desse porte é cara em memória e muito cara para eliminação de Gauss.

O artigo usa GMRES porque o sistema aumentado com multiplicadores de Lagrange é grande e não é positivo definido.

### 5.2. Estrutura esparsa mínima

Implementar primeiro formato COO didático:

```c
typedef struct {
    int row;
    int col;
    double value;
} SparseEntry;
```

Depois converter para CSR:

```text
row_ptr
col_ind
values
```

### 5.3. Tarefas

- [x] Criar módulo `sparse_matrix`.
- [x] Implementar estrutura COO.
- [x] Permitir acumular entradas repetidas.
- [x] Implementar ordenação por `(row, col)`.
- [x] Implementar compressão de duplicatas.
- [x] Implementar conversão COO → CSR.
- [x] Implementar produto matriz-vetor CSR.
- [x] Criar testes unitários para:

  - inserção;
  - duplicatas;
  - compressão;
  - produto matriz-vetor;
  - comparação com matriz densa pequena.

---

## 6. Prioridade 4 — Montagem esparsa da rigidez

### 6.1. Objetivo

Criar equivalente esparso de:

```c
global_stiffness_assemble_dense(...)
```

sem remover a versão densa.

### 6.2. Tarefas

- [x] Criar `global_stiffness_assemble_sparse_coo(...)`.
- [x] Reutilizar a montagem local já existente.
- [x] Acumular entradas locais em COO.
- [x] Comprimir duplicatas.
- [x] Converter para CSR.
- [x] Comparar, em problemas pequenos, a matriz esparsa com a matriz densa.
- [x] Criar teste:

  ```text
  sparse K ≈ dense K
  ```

para `3x3x3`.

---

## 7. Prioridade 5 — Sistema aumentado esparso

### 7.1. Objetivo

Montar em formato esparso:

```text
[ K  G^T ]
[ G   0  ]
```

### 7.2. Tarefas

- [x] Criar montagem esparsa de `G`.
- [x] Criar montagem esparsa do sistema aumentado.
- [x] Comparar com o sistema denso em caso pequeno.
- [x] Criar teste:

  ```text
  sparse A_aug * x ≈ dense A_aug * x
  ```

---

## 8. Prioridade 6 — GMRES

### 8.1. Objetivo

Implementar GMRES para resolver:

```text
A_aug x = b_aug
```

com matriz esparsa CSR.

### 8.2. Tarefas

- [x] Implementar GMRES básico sem precondicionador.
- [x] Criar teste com sistema pequeno conhecido.
- [x] Testar contra `dense_solve(...)` em sistema pequeno.
- [x] Registrar:

  - número de iterações;
  - resíduo inicial;
  - resíduo final;
  - tolerância;
  - convergência ou falha.
- [x] Implementar versão reiniciada, por exemplo GMRES(m), se necessário.
- [ ] Avaliar necessidade de precondicionador.

### 8.3. Possíveis precondicionadores futuros

- [ ] Jacobi.
- [ ] ILU simples.
- [ ] Escalonamento diagonal.

---

## 9. Prioridade 7 — Reprodução do cubo com refinamento suficiente

### 9.1. Objetivo

Resolver o problema do cubo com refinamento suficiente para aproximar a figura do artigo.

### 9.2. Casos-alvo

Casos uniformes iniciais:

```text
11x11x11 nós, 15x15x15 células
13x13x13 nós, 15x15x15 células
15x15x15 nós, 15x15x15 células
```

Depois, testar nuvens não uniformes somente após diagnóstico de conectividade e condicionamento.

### 9.3. Métricas

- [ ] erro relativo discreto global;
- [ ] erro relativo discreto interno;
- [ ] erro máximo absoluto;
- [ ] erro no plano `x = 5.33`;
- [ ] erro em regiões próximas às arestas superiores;
- [ ] número de iterações GMRES;
- [ ] tempo de montagem;
- [ ] tempo de solução;
- [ ] memória aproximada.

---

## 10. Prioridade 8 — Reprodução da figura

### 10.1. Figura de interesse

A figura principal deve ser reproduzida em plano semelhante ao usado no artigo, especialmente:

```text
x = 5.33
```

### 10.2. Tarefas

- [ ] Criar app para exportar plano `x = 5.33`.
- [ ] Gerar CSV:

  ```text
  y,z,V_num,V_exact,abs_error
  ```

- [ ] Criar script Python para plotar:

  - solução numérica;
  - solução analítica;
  - erro absoluto;
  - mapa de contorno.

- [ ] Comparar erro máximo com valor reportado no artigo.
- [ ] Documentar divergências.

---

## 11. Prioridade 9 — Nuvens não uniformes

### 11.1. Status

A primeira tentativa de nuvem não uniforme foi útil como experimento, mas não deve ser a linha principal agora.

### 11.2. Manter como investigação secundária

- [ ] Reavaliar depois da implementação esparsa.
- [ ] Reavaliar depois do diagnóstico de conectividade.
- [ ] Reavaliar depois do diagnóstico de condicionamento MLS.
- [ ] Evitar nuvens dominadas por pontos de superfície.
- [ ] Buscar equilíbrio entre:

  - nós internos;
  - nós de superfície;
  - refinamento próximo às arestas superiores;
  - estabilidade da matriz de momento.

---

## 12. O que estava sendo proposto antes, mas deve ficar pausado

As seguintes ideias estavam no caminho, mas devem ser pausadas:

- [ ] Comparar estatisticamente CSV regular versus não uniforme.
- [ ] Ajustar empiricamente a nuvem não uniforme.
- [ ] Criar muitos gráficos de erro por região.
- [ ] Refinar heurísticas de concentração de nós no topo.
- [ ] Investigar overshoot em nuvens não uniformes.

Essas tarefas são úteis, mas são secundárias. O gargalo principal agora é estrutural:

```text
conectividade MLS
condicionamento
matriz esparsa
GMRES
refinamento suficiente
```

---

## 13. Ordem recomendada de execução a partir de agora

1. Diagnóstico de conectividade MLS.
2. Diagnóstico da matriz de momento MLS.
3. Matriz esparsa COO/CSR.
4. Montagem esparsa de K.
5. Montagem esparsa de G.
6. Sistema aumentado esparso.
7. GMRES básico.
8. Comparação GMRES versus solver denso em casos pequenos.
9. Rodar casos uniformes maiores.
10. Exportar plano `x = 5.33`.
11. Reproduzir figura do artigo.
12. Só depois voltar para nuvens não uniformes.

---

## 14. Critério de sucesso intermediário

Antes de buscar a figura final, o projeto deve conseguir rodar:

```text
15x15x15 nós
15x15x15 células de integração
Dirichlet via Lagrange
GMRES
exportação de plano x = 5.33
comparação com solução analítica
```

com relatório contendo:

```text
nodes
constraints
nonzeros
GMRES iterations
residual
relative error
max abs error
runtime
```

---

## 15. Critério de sucesso final

O projeto será considerado cientificamente consistente quando conseguir:

- reproduzir qualitativamente a distribuição de potencial no cubo;
- obter erro interno decrescente com refinamento;
- exportar figura comparável à do artigo;
- documentar diferenças em relação ao artigo;
- explicar claramente o papel de:

  - suporte nodal;
  - cardinalidade local;
  - matriz de momento MLS;
  - multiplicadores de Lagrange;
  - GMRES;
  - distribuição regular versus não uniforme de nós.
