# GEMINI_REVIEW — Revisão local, didática e documental

## Data

2026-05-06

## Objetivo da revisão

Revisar a implementação recente do `reproduce_cube_sparse` feita pelo Codex, com foco em:

- clareza;
- manutenção;
- mensagens de erro;
- comentários didáticos;
- documentação;
- consistência com o `TODO.md` e a auditoria do Claude.

## Arquivos revisados

- [x] `TODO.md`
- [x] `docs/dev/AI_HANDOFF.md`
- [x] `docs/dev/CLAUDE_AUDIT.md`
- [x] `docs/dev/CODEX_IMPLEMENTATION_LOG.md`
- [x] `docs/dev/TEST_REPORT.md`
- [x] `apps/reproduce_cube_sparse.c`
- [x] `CMakeLists.txt`
- Nenhum arquivo em `include/` ou `src/` foi modificado na rodada anterior, então a revisão focou no app.

## Pontos verificados

### 1. Clareza do app (`reproduce_cube_sparse.c`)

- [x] **Parâmetros físicos e numéricos claros:** `L`, `V0`, contagens de nós/células e parâmetros do GMRES são bem definidos e passados para a função `run_case`.
- [x] **Relatório terminal compreensível:** A saída do console é bem estruturada, separando diagnóstico, montagem, solução e erro. O bloco `--- Required report ---` é consistente com as métricas obrigatórias.
- [x] **Tratamento de erro adequado:** O app verifica falhas de alocação e usa `goto` para limpeza, um padrão C robusto. Os critérios de parada (`support_lt_4 > 0`, `mls_failures > 0`, `gmres_converged = NO`) estão implementados e funcionam como "disjuntores" de segurança, conforme exigido.
- [x] **Separação entre caso pequeno e caso principal:** A CLI com `--case sanity` e `--case target` e a lógica na função `main` separam claramente o caso de validação (5x5x5 com comparação densa) do caso de produção (11x11x11, GMRES apenas).

Observações:

```text
O código do app é de alta qualidade, robusto e fácil de seguir. As decisões de implementação do Codex estão bem alinhadas com as auditorias e os requisitos do projeto.
```

### 2. Comentários didáticos

Foram adicionados comentários para explicar o "porquê" de certas decisões numéricas, melhorando o valor didático do código.

- [x] **Base MLS e cardinalidade:** Adicionado comentário explicando por que a cardinalidade mínima é 4 para a base linear `[1, x, y, z]` e por que menos de 8 nós é uma região suspeita.
- [x] **Sistema aumentado `[K G^T; G 0]`:** Adicionado comentário explicando que essa estrutura vem dos multiplicadores de Lagrange para impor as condições de Dirichlet.
- [x] **Uso de GMRES:** Adicionado comentário explicando que GMRES é necessário porque o sistema aumentado é simétrico, mas indefinido, tornando o Gradiente Conjugado inadequado.
- [x] **Erro máximo:** Adicionado comentário na seção de erros para explicar por que um `max_abs_error` da ordem de `V0` é esperado e não é um bug, destacando a importância do erro relativo interior.

Observações:

```text
Os novos comentários conectam diretamente o código às discussões teóricas presentes nos documentos de design (`CLAUDE_AUDIT.md`, `TODO.md`), cumprindo um dos objetivos principais do projeto.
```

### 3. Documentação

- [x] **`TODO.md` atualizado:** O `TODO.md` reflete com precisão o estado do projeto. As tarefas de diagnóstico, esparsificação e GMRES estão marcadas como concluídas, e a próxima prioridade (refinamento) está corretamente identificada.
- [x] **`TEST_REPORT.md` e `CODEX_IMPLEMENTATION_LOG.md` suficientes:** Ambos os relatórios são extremamente detalhados e fornecem um registro claro do que foi feito, como foi testado e quais foram os resultados. O problema de ambiente com `ctest` foi bem documentado.
- [x] **Nuvens não uniformes continuam pausadas:** Todos os documentos relevantes (`TODO.md`, `CLAUDE_AUDIT.md`) confirmam que esta linha de trabalho está corretamente em pausa.
- [x] **Próximo passo está claro:** A recomendação é avançar para malhas mais refinadas (13x13x13, 15x15x15) usando o pipeline esparso validado.

Observações:

```text
A documentação do projeto é exemplar, permitindo uma rastreabilidade clara entre requisitos, auditorias, implementação e resultados.
```

## Pequenas correções feitas

- Em `apps/reproduce_cube_sparse.c`, o campo `relative_error_internal` foi renomeado para `relative_error_interior` na struct `CubeSparseReport` e nas funções de relatório para maior consistência com a nomenclatura usada na saída do console.
- Adicionados os comentários didáticos mencionados acima.

## Problemas encontrados

- Nenhum problema bloqueador foi encontrado.
- A única inconsistência foi o nome do campo de erro (`internal` vs. `interior`), que foi corrigido.
- O problema com o `ctest` é específico do ambiente de execução e não um bug do projeto, mas foi bem documentado no `TEST_REPORT.md`.

## Riscos restantes

Os riscos são os mesmos identificados na auditoria do Claude e permanecem gerenciáveis para os próximos passos:

1.  **Memória (Alto para >17x17x17):** A montagem ainda passa por uma matriz `K` densa intermediária. Isso é aceitável para 13x13x13 e 15x15x15, mas se tornará um gargalo para refinamentos muito maiores.
2.  **Convergência do GMRES (Baixo):** Os parâmetros de `restart` e `max_iter` podem precisar de ajuste para os casos maiores, conforme recomendado pela auditoria.
3.  **Condicionamento do MLS (Baixo):** O condicionamento da matriz de momento `A(x)` pode piorar com o refinamento. Isso deve ser monitorado através do campo `max_cond_estimate` nos relatórios.

## Recomendação

```text
Pode avançar para 13x13x13 e 15x15x15.
```

Justificativa:

```text
O pipeline esparso está estável, validado e robusto. Os testes para 5x5x5 e 11x11x11 passaram com sucesso, e os critérios de parada (conectividade, falhas de MLS, convergência do GMRES) foram satisfeitos. Os riscos são conhecidos e monitoráveis.

A recomendação é aplicar as pequenas melhorias de clareza e comentários propostas nesta revisão e, em seguida, prosseguir com a execução dos casos de refinamento `13x13x13` e `15x15x15`, ajustando os parâmetros do GMRES conforme a auditoria do Claude sugeriu. Não há necessidade de revisões adicionais de arquitetura ou correção de bugs antes disso.
```