# efg-3d-electromagnetics-2006

Reprodução didática em C11 do artigo de Parreira et al. (2006) sobre o método Element-Free Galerkin aplicado a problemas eletromagnéticos tridimensionais, com tradução comentada, formulação matemática, testes numéricos e comparação com FEM.

## Estado atual

O projeto ainda está na fase de esqueleto. A estrutura em C11 já compila, mas o solver EFG completo ainda não foi implementado.

Esta primeira etapa contém apenas módulos vazios e didáticos para organizar o trabalho futuro:

- `vec3`
- `weight`
- `gauss`
- `analytical`
- `error_norm`
- `nodes`
- `mls`
- `assembly`
- `gmres`

Os arquivos de código incluem TODOs apontando para a tradução em `docs/`, para evitar a introdução de fórmulas fora das referências do artigo.

## Estrutura

- `include/efg/`: headers públicos dos módulos.
- `src/`: fontes dos módulos.
- `apps/`: executáveis simples, começando por `smoke_test`.
- `tests/`: espaço reservado para testes automatizados.
- `data/input/`: entradas futuras.
- `data/output/`: saídas numéricas futuras em CSV.
- `scripts/`: scripts futuros de validação e apoio.

## Como compilar

```sh
cmake -S . -B build
cmake --build build
ctest --test-dir build
```

O executável `smoke_test` apenas verifica se todos os headers e fontes do esqueleto compilam e linkam corretamente. Ele não executa o método EFG.

## Navegação

[Capa e resumo](docs/00_capa_resumo.md) |
[1. Introdução](docs/01_introducao.md) |
[2. Formulação matemática](docs/02_formulacao_matematica.md) |
[3. Resultados numéricos](docs/03_resultados_numericos.md) |
[4. Conclusão e referências](docs/04_conclusao_referencias.md)

> **Fonte:** tradução técnica do artigo *The Element-Free Galerkin Method in Three-Dimensional Electromagnetic Problems*, publicado em *IEEE Transactions on Magnetics*, v. 42, n. 4, abril de 2006. DOI: [10.1109/TMAG.2006.872014](https://doi.org/10.1109/TMAG.2006.872014).
