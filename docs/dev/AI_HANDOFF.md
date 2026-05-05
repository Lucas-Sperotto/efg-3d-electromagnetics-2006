# AI_HANDOFF — Cooperação entre Claude Code, Codex e Gemini Code Assist

Este arquivo organiza a cooperação entre IAs no desenvolvimento do projeto
`efg-3d-electromagnetics-2006`.

O objetivo do projeto é reproduzir, de forma didática e verificável, o problema
eletrostático 3D apresentado no artigo:

> The Element-Free Galerkin Method in Three-Dimensional Electromagnetic Problems  
> Parreira, Silva, Fonseca e Mesquita, IEEE Transactions on Magnetics, 2006.

## Papéis das IAs

### Claude Code

Papel principal:

- auditor matemático;
- arquiteto técnico;
- revisor de coerência científica;
- verificador de riscos antes de avançar para casos maiores.

Claude deve priorizar:

- formulação EFG/MLS;
- cardinalidade mínima da base;
- suporte nodal;
- condicionamento da matriz de momento MLS;
- montagem do sistema aumentado;
- coerência entre implementação e artigo.

Claude não deve priorizar:

- grandes refatorações;
- otimizações prematuras;
- mudanças estéticas extensas.

---

### Codex

Papel principal:

- implementação em C/CMake;
- criação de apps;
- criação de testes;
- execução de build e ctest;
- registro objetivo dos resultados.

Codex deve priorizar:

- alterações pequenas e testáveis;
- preservação da versão densa como referência;
- implementação esparsa;
- GMRES;
- relatórios de execução.

Codex não deve priorizar:

- documentação longa;
- reescrita ampla da arquitetura;
- mudança em nuvens não uniformes sem autorização.

---

### Gemini Code Assist

Papel principal:

- revisão local no VS Code;
- clareza de código;
- comentários didáticos;
- organização de documentação;
- detecção de inconsistências simples.

Gemini deve priorizar:

- legibilidade;
- nomes de variáveis;
- mensagens de erro;
- documentação dos apps;
- coerência do TODO.

Gemini não deve priorizar:

- grandes mudanças matemáticas;
- alteração de APIs públicas;
- refatoração profunda sem necessidade.

---

## Estado técnico atual

Já existem ou estão em desenvolvimento:

- diagnóstico MLS;
- estrutura esparsa COO/CSR;
- montagem esparsa da matriz de rigidez;
- montagem esparsa do sistema aumentado com multiplicadores de Lagrange;
- GMRES reiniciado;
- aplicação para reprodução do cubo eletrostático.

## Prioridade atual

A prioridade atual é validar o pipeline esparso no problema regular do cubo:

- cubo com `L = 10`;
- potencial superior `V0 = 10`;
- grade regular inicial `11x11x11`;
- integração `15x15x15`;
- quadratura Gauss `2x2x2`;
- Dirichlet via multiplicadores de Lagrange;
- solução com GMRES.

## Critérios de bloqueio

Não avançar para casos maiores se ocorrer:

```text
support_lt_4 > 0
mls_failures > 0
gmres_converged = false
residual_final alto sem explicação
erro interno incompatível com a versão densa
```

## Métricas obrigatórias nos relatórios

Sempre que o app principal for executado, registrar:

```text
nodes
constraints
augmented_size
K_nnz
A_aug_nnz
support_min
support_mean
support_max
support_lt_4
support_lt_8
mls_failures
gmres_converged
gmres_iterations
residual_initial
residual_final
relative_error_global
relative_error_internal
max_abs_error
assembly_time_s
solve_time_s
```

## O que está pausado

Até validação do pipeline esparso, manter pausado:

* nuvens não uniformes;
* octree;
* precondicionadores;
* novos gráficos;
* ajuste empírico de distribuição de nós;
* reprodução estética final da figura.

## Fluxo recomendado

1. Claude audita.
2. Codex implementa.
3. Gemini revisa.
4. Claude valida se pode avançar.
