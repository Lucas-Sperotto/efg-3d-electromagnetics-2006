# 3. Resultados numéricos

Para avaliar a correção e a precisão do método, foi utilizado um problema com solução analítica viável.

## 3.1. Problema de referência: caixa cúbica eletrostática

Considere uma caixa cúbica com paredes condutoras perfeitas, como mostrado na Fig. 1. Para caracterizar esse problema como eletromagnético, as paredes condutoras com diferentes potenciais não podem estar conectadas. Ele pode ser modelado como uma caixa com um isolante infinitesimal entre as paredes laterais e a parede superior.

**Fig. 1.** Caixa cúbica eletrostática.

Esse problema é um bom *benchmark* porque possui solução analítica e porque o potencial eletrostático apresenta fortes variações nos cantos superiores. A solução analítica para $L = 10$ e $V_0 = 10\,\text{V}$ é

$$
V = \sum_{\substack{n\ \text{ímpar}}} \sum_{\substack{m\ \text{ímpar}}}
a_{mn}\,\sin\left(\frac{n\pi}{10}x\right)\,
\sin\left(\frac{m\pi}{10}y\right)\,
\sinh(k_{mn}z),
$$

em que

$$
a_{mn} = \frac{160}{mn\pi^2\sinh(10k_{mn})}
$$

e

$$
k_{mn} = \sqrt{\left(\frac{n\pi}{10}\right)^2 + \left(\frac{m\pi}{10}\right)^2}.
$$

## 3.2. Distribuição de nós e domínio de influência

Para resolver esse problema com o EFG, foi utilizada uma distribuição não uniforme de nós, com nós densamente alocados nos cantos superiores da caixa e distribuídos de forma mais grosseira ao longo da região em que os potenciais apresentam comportamento mais regular. Essa distribuição pode ser vista na Fig. 2.

**Fig. 2.** Distribuição dos nós para o problema de referência.

O tamanho do domínio de influência $d_{mI}$ deve ser definido com cuidado. Se esse tamanho for muito grande, muitos nós serão incluídos no domínio de influência, de modo que o custo computacional da construção das funções de forma aumentará de maneira proibitiva. Por outro lado, se $d_{mI}$ for muito pequeno, o domínio poderá não ser completamente coberto e as funções de forma não serão construídas corretamente.

Assim, um bom tamanho para o domínio de influência depende da densidade local de nós em cada região. Por isso, foi utilizada uma técnica de pré-processamento que percorre os nós e realiza uma busca pelos $k$ vizinhos mais próximos (*k-nearest neighbor*) em uma árvore que contém esses nós.

Com a informação dos nós mais próximos, é possível estimar quão densa é a distribuição de nós naquela localização. Então, o tamanho do domínio de influência — e, consequentemente, o suporte da função de forma — pode ser corretamente definido para cada nó. Essa busca pode ser realizada de maneira eficiente e não exige tempo de processamento considerável.

Outro aspecto importante, especialmente em dimensões superiores, é que, à medida que o domínio de influência diminui, a esparsidade da matriz de rigidez aumenta, pois o número de sobreposições entre funções de forma é reduzido. Isso conduz a uma solução mais rápida do sistema linear e a um menor uso de memória.

Neste trabalho, foram localizados os três vizinhos mais próximos de cada nó, e $d_{mI}$ foi definido como a distância do nó até o terceiro vizinho multiplicada por $1{,}5$. Assim, garante-se que pelo menos quatro funções de forma cubram cada região do domínio.

## 3.3. Integração numérica e comparação com FEM

A integração foi realizada utilizando $15 \times 15 \times 15$ células de integração. A quadratura de Gauss de ordem 2 foi suficiente para obter uma boa aproximação, apesar da complexidade das funções de forma. Foi adotada a função peso spline cúbica da Eq. (10), com construção cúbica do domínio de influência.

Ao utilizar a função delta de Dirac como função interpoladora dos multiplicadores de Lagrange, obtém-se uma simplificação, pois duas integrais de superfície da Eq. (7) não precisam ser realizadas.

As linhas de contorno da solução foram plotadas no plano $yz$ em $x = 5{,}33$, como mostrado na Fig. 3. Esse problema também foi resolvido usando FEM com os mesmos nós em uma malha de Delaunay. A norma do erro relativo dos valores do potencial, calculada pela Eq. (16), foi de $3{,}56\%$ no FEM e de $1{,}40\%$ no EFG. O erro máximo encontrado foi de $1{,}08\,\text{V}$ no FEM e de $0{,}15\,\text{V}$ no EFG.

**Fig. 3.** Gráficos de contorno do potencial em $x = 5{,}33$.

A norma relativa do erro é dada por

$$
e = \left[
\frac{\displaystyle\int_{\Omega}\left(V_{\text{EFG}} - V_{\text{exato}}\right)^2\,d\Omega}
{\displaystyle\int_{\Omega}V_{\text{exato}}^2\,d\Omega}
\right]^{1/2}.
$$

A comparação entre as soluções por EFG e FEM não considerou os tempos computacionais. É evidente que a geração das matrizes do sistema é mais custosa para o EFG do que para o FEM. Entretanto, os tempos de pré-processamento e pós-processamento são significativamente reduzidos, uma vez que nenhuma malha é necessária.

## 3.4. Convergência com o refinamento dos nós

Para analisar a convergência do método, foi utilizada uma distribuição uniforme de nós. A curva de convergência do método pode ser observada na Fig. 4, em que $d$ indica a distância entre nós vizinhos.

**Fig. 4.** Norma relativa do erro durante o refinamento da distribuição de nós.

Como pode ser observado, foi obtida uma taxa de convergência de primeira ordem, pois as funções de forma escolhidas nas Eqs. (8) e (9) são interpolantes de primeira ordem.

A solução obtida com duas funções peso foi avaliada: a função peso spline cúbica da Eq. (10) e a função peso quadrática da Eq. (11). Embora a primeira seja mais complexa do que a segunda, o erro gerado com a função quadrática foi ligeiramente menor do que aquele obtido com a spline cúbica. Entretanto, a diferença observada no erro entre as duas funções peso indica que o uso de qualquer uma delas não deve representar uma dificuldade.

## 3.5. Influência do número de células de integração

A Fig. 5 mostra o comportamento da solução em relação à integração numérica. Ao utilizar o método de Gauss-Legendre, a grade de integração é criada dividindo todo o domínio em caixas cúbicas iguais, denominadas células de Gauss.

**Fig. 5.** Norma relativa do erro ao aumentar o número de células de integração.

A integração das derivadas das funções de forma ainda é um problema em aberto no EFG, pois elas não possuem forma analítica. Ao aumentar o número de células, observa-se que a solução converge rapidamente para a solução analítica, e aumentos adicionais não conduzem a melhorias significativas.

Com base nesses resultados, um parâmetro pode ser definido para ajustar a densidade das células de integração à densidade dos nós. Trabalhando de acordo com essa estratégia, tempo computacional não será gasto onde ele não é necessário.

Acoplar o número de células de integração à densidade de nós é uma estratégia que pode ser prontamente adotada de duas maneiras. Uma possibilidade é utilizar uma estrutura de dados, como uma *octree*, para definir regiões em que a integração deve ser mais fina ou mais grosseira de acordo com a densidade dos nós, estabelecendo então o número de células para cada região. Como segunda opção, a ordem de integração pode ser alterada de célula para célula, de modo a levar em conta a densidade de nós no interior das células.
