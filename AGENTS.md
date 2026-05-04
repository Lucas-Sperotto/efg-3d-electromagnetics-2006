# AGENTS.md

Este repositório reproduz didaticamente, em linguagem C, o artigo
"The Element-Free Galerkin Method in Three-Dimensional Electromagnetic Problems".

## Regras gerais

1. O código deve ser escrito em C, preferencialmente C11.
2. O objetivo principal é didático, não otimização prematura.
3. Toda função que implemente uma equação deve citar no comentário:
   - número da equação;
   - arquivo `.md` de referência;
   - significado físico/matemático.
4. Não alterar fórmulas sem consultar os arquivos da pasta `docs/`.
5. Toda saída numérica deve ser salva em CSV.
6. Não declarar uma figura como reproduzida sem comparação com:
   - solução analítica;
   - dados do artigo;
   - ou gráfico correspondente.
7. Preferir funções pequenas, testáveis e bem comentadas.
8. Separar claramente:
   - código matemático básico;
   - montagem do método EFG;
   - resolução do sistema linear;
   - scripts de validação.
9. Não usar bibliotecas externas na primeira versão, exceto a biblioteca padrão de C e `math.h`.
10. Futuramente, bibliotecas externas podem ser introduzidas para matrizes esparsas ou GMRES, mas somente depois da versão didática.
