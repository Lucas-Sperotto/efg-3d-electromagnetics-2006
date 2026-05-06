# GEMINI_REVIEW — Revisão local e documental

## Data

2026-05-08

## Objetivo da revisão

Revisar a implementação do script de plotagem `scripts/plot_cube_plane.py` e a consistência dos resultados gerados, validando que o projeto está pronto para a comparação final com a Figura 3 do artigo.

## Arquivos revisados

- [x] `TODO.md`
- [x] `docs/dev/AI_HANDOFF.md`
- [x] `docs/dev/CLAUDE_AUDIT.md`
- [x] `docs/dev/CODEX_IMPLEMENTATION_LOG.md`
- [x] `docs/dev/TEST_REPORT.md`
- [x] `apps/reproduce_cube_sparse.c`
- [x] `scripts/plot_cube_plane.py`
- [x] `CMakeLists.txt`

## Pontos verificados

### 1. Pipeline Numérico (`reproduce_cube_sparse.c`)

- [x] **Refinamento concluído:** Os casos `refine13` e `refine15` foram executados com sucesso, mostrando a convergência esperada do erro relativo interior.
- [x] **Exportação do plano concluída:** O caso `--case plane15` executa a simulação 15x15x15 e exporta corretamente o arquivo `data/output/cube_plane_x_5_33_refine15.csv` com os dados necessários para a plotagem.

### 2. Script de Plotagem (`scripts/plot_cube_plane.py`)

- [x] **Funcionalidade:** O script lê o CSV, processa os dados e gera os quatro gráficos de contorno solicitados (`V_num`, `V_exact`, `abs_error`, `comparison`).
- [x] **Robustez:** O script valida o CSV de entrada (colunas, valores finitos, duplicatas) e trata erros de forma clara.
- [x] **Métricas:** Além dos gráficos, o script calcula e salva métricas importantes, incluindo o erro relativo no **interior** do plano, que é a métrica mais relevante para a análise de convergência, isolando os efeitos de borda.
- [x] **Saídas:** As saídas (`.png`, `.txt`, `.csv`) são salvas no diretório `data/output/figures/`, mantendo a organização do projeto.

### 3. Documentação

- [x] **`TODO.md`:** O arquivo reflete o estado atual. As tarefas de refinamento, exportação e plotagem estão marcadas como concluídas. A próxima prioridade é a comparação com o artigo.
- [x] **`TEST_REPORT.md`:** O relatório documenta a execução bem-sucedida do caso `plane15` e a geração do CSV, validando a entrada para o script de plotagem.

## Problemas encontrados

- Nenhum problema bloqueador foi encontrado.
- O pipeline completo, da simulação à visualização, está funcional e robusto.

## Riscos restantes

Os riscos são os mesmos identificados nas auditorias anteriores e permanecem gerenciáveis, mas não impactam a próxima etapa de análise:

1.  **Memória (Alto para >17x17x17):** A montagem via matriz densa intermediária continua sendo um gargalo para refinamentos muito maiores, mas o alvo de 15x15x15 foi atingido com sucesso.
2.  **Condicionamento do MLS (Baixo):** O condicionamento pode piorar com refinamentos futuros, mas se manteve estável e dentro dos limites aceitáveis (`max_cond_estimate < 1e6`) para os casos executados.

## Recomendação

```text
O projeto está pronto para a etapa final de análise científica. A recomendação é avançar para a comparação visual dos gráficos gerados com a Figura 3 do artigo e documentar as conclusões.
```

Justificativa:

```text
O pipeline esparso está estável, validado e robusto. Os testes para 5x5x5 e 11x11x11 passaram com sucesso, e os critérios de parada (conectividade, falhas de MLS, convergência do GMRES) foram satisfeitos. Os riscos são conhecidos e monitoráveis.

A recomendação é aplicar as pequenas melhorias de clareza e comentários propostas nesta revisão e, em seguida, prosseguir com a execução dos casos de refinamento `13x13x13` e `15x15x15`, ajustando os parâmetros do GMRES conforme a auditoria do Claude sugeriu. Não há necessidade de revisões adicionais de arquitetura ou correção de bugs antes disso.
```