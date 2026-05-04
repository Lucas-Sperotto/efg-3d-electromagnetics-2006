# Condições de contorno e multiplicadores de Lagrange

Este documento formaliza a baseline conceitual para a implementação futura das condições de contorno do problema do cubo e da imposição de Dirichlet no método EFG.

Ele não introduz solver EFG nem novas fórmulas além das relações já presentes na tradução e das simplificações explicitamente indicadas como interpretação de implementação.

---

## Condições de contorno do cubo

Para o problema de referência da caixa cúbica eletrostática, considera-se um cubo de lado `L` com potencial prescrito na face superior e potencial nulo nas demais faces.

As condições de contorno do benchmark ficam:

- `V = 0` em `x = 0`
- `V = 0` em `x = L`
- `V = 0` em `y = 0`
- `V = 0` em `y = L`
- `V = 0` em `z = 0`
- `V = V0` em `z = L`

Com os valores numéricos do artigo:

- `L = 10`
- `V0 = 10 V`

TODO: confirmar no artigo original se há alguma convenção adicional para o isolante infinitesimal entre a parede superior e as paredes laterais.

---

## Por que Dirichlet não é imposto diretamente no EFG

No FEM nodal usual, as funções de forma possuem a propriedade delta de Kronecker. Isso permite associar diretamente um grau de liberdade nodal ao valor da solução naquele nó.

No EFG, conforme `docs/02_formulacao_matematica.md`, as funções de forma não satisfazem essa propriedade. Assim, o coeficiente associado a um nó não é simplesmente o valor físico do potencial naquele ponto.

Consequência prática:

- Não se deve impor `V = Vd` alterando diretamente graus de liberdade nodais como no FEM nodal clássico.
- A condição essencial de Dirichlet precisa de uma técnica adicional de imposição.
- O artigo escolhe multiplicadores de Lagrange para evitar parâmetros empíricos, como os que aparecem no método da penalidade.

---

## Ideia dos multiplicadores de Lagrange

A tradução em `docs/02_formulacao_matematica.md`, seção 2.1, introduz um funcional energético modificado com multiplicadores de Lagrange.

A ideia é acrescentar uma restrição ao problema variacional para que a solução satisfaça a condição de Dirichlet na fronteira prescrita.

Na forma discreta, isso produz um sistema em blocos com:

- matriz de rigidez `K`;
- bloco de acoplamento de fronteira `G`;
- transposta `G^T`;
- incógnitas do potencial `V`;
- incógnitas dos multiplicadores `lambda`;
- vetor volumétrico `f`;
- vetor de contorno `b`.

Esse sistema deixa de ser definido positivo, como a própria tradução observa, e por isso o artigo usa GMRES para resolvê-lo.

TODO: revisar a dimensão exata de `lambda`, `N_lambda` e `N_i` antes da implementação.

---

## Simplificação com delta de Dirac para os multiplicadores

Em `docs/03_resultados_numericos.md`, seção 3.3, a tradução afirma que, ao utilizar a função delta de Dirac como função interpolante dos multiplicadores de Lagrange, duas integrais de superfície não precisam ser realizadas.

Interpretação para implementação, a ser revisada manualmente:

- o termo `G_ij` passa a ser avaliado pontualmente nos pontos de fronteira de Dirichlet;
- o termo `b_i` passa a ser o valor prescrito no ponto de fronteira correspondente.

Essa simplificação deve ser tratada com cuidado porque a documentação atual ainda não explicita a expressão resultante nem a numeração exata da equação referenciada.

TODO: confirmar contra o artigo original se a avaliação pontual deve usar exatamente o nó de fronteira, ponto de quadratura de fronteira ou outro conjunto de pontos.

TODO: definir no código uma estrutura explícita para pontos de Dirichlet antes de montar `G` e `b`.

---

## Regras para a implementação futura

- Implementar a fronteira do cubo antes de montar os termos de Lagrange.
- Separar claramente dados geométricos, valores prescritos e montagem algébrica.
- Não misturar a simplificação por delta de Dirac com uma integração de superfície genérica na mesma função.
- Documentar em cada função C o arquivo Markdown de referência e a equação correspondente.
- Marcar como TODO qualquer detalhe que ainda dependa de conferência contra o artigo original.
