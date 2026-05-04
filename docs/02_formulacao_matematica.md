# 2. Formulação matemática

## 2.1. Forma fraca e sistema linear

Problemas eletrostáticos podem ser descritos pela forma forte definida pela equação de Poisson:

$$
\vec{\nabla} \cdot \left(\varepsilon \vec{\nabla} V\right) = -\rho
$$

com as condições de contorno:

$$
V = V_d \quad \text{em } \Gamma_d
\qquad \text{e} \qquad
-\varepsilon \frac{\partial V}{\partial n} = 0 \quad \text{em } \Gamma_n .
$$

Em que $V$ é o potencial eletrostático, $\varepsilon$ é a permissividade elétrica, $\rho$ é a densidade volumétrica de carga, $\Gamma_n$ é a fronteira com condição de Neumann homogênea e $V_d$ é o potencial imposto na fronteira de Dirichlet $\Gamma_d$.

As funções de forma do EFG não satisfazem a propriedade do delta de Kronecker. Portanto, as condições de contorno essenciais não podem ser impostas diretamente, como é feito no FEM e no método das diferenças finitas. Uma técnica separada deve ser utilizada para impô-las.

Entre as técnicas mais utilizadas, podem ser mencionados o método dos multiplicadores de Lagrange, o método da penalidade, o acoplamento com elementos finitos e o método de Nitsche [6]. Neste trabalho, foi escolhido o método dos multiplicadores de Lagrange, pois não exige a definição de parâmetros adicionais, às vezes empíricos, como ocorre no método da penalidade, nem introduz características híbridas, como a conectividade de elementos, presente no acoplamento FEM-EFG. Além disso, sua implementação é direta e bons resultados são obtidos com seu uso.

Para utilizar os multiplicadores de Lagrange, emprega-se um funcional energético modificado, como mostrado a seguir:

$$
\Pi(V,\lambda) =
\frac{1}{2}\int_{\Omega} \varepsilon \vec{\nabla} V \cdot \vec{\nabla} V d\Omega -
\int_{\Omega} \rho V d\Omega +
\int_{\Gamma_d} \lambda (V - V_d)\,d\Gamma .
$$

A avaliação do ponto estacionário do funcional em relação às variações em $V$ e em $\lambda$ gera a seguinte forma fraca:

**Encontrar $V \in H^1(\Omega)$ e $\lambda \in H^{-1/2}(\Gamma_d)$ tais que:**

$$
\int_{\Omega} \varepsilon \vec{\nabla} V \cdot \vec{\nabla} u d\Omega +
\int_{\Gamma_d} \lambda u d\Gamma =
\int_{\Omega} \rho u d\Omega,
\qquad
\forall u \in H^1(\Omega),
$$

$$
\int_{\Gamma_d} \gamma (V - V_d)\,d\Gamma = 0,
\qquad
\forall \gamma \in H^{-1/2}(\Gamma_d),
$$

em que $u$ e $\gamma$ são funções teste.

O método de Galerkin é utilizado para obter as equações discretas da forma fraca por meio das seguintes aproximações:

$$
u(x) \approx \sum_{i=1}^{N} \Phi_i(x)u_i,
$$

$$
V(x) \approx \sum_{i=1}^{N} \Phi_i(x)V_i,
$$

$$
\lambda(x) \approx \sum_{i=1}^{N_{\lambda}} N_i(x)\lambda_i .
$$

Em que $x=(x,y,z)$, $\Phi_i(x)$ são as funções de forma do EFG e $N_i(x)$ são os interpolantes dos multiplicadores de Lagrange. O seguinte sistema linear é obtido:

$$
\begin{bmatrix}
K & G^T \\
G & 0
\end{bmatrix}
\begin{bmatrix}
V \\
\lambda
\end{bmatrix} =
\begin{bmatrix}
f \\
b
\end{bmatrix} .
$$

Com:

$$
K_{ij}=\int_{\Omega} \varepsilon \nabla \Phi_i \cdot \nabla \Phi_j\,d\Omega,
$$

$$
G_{ij}=\int_{\Gamma_d} N_i\Phi_j\,d\Gamma,
$$

$$
b_i=\int_{\Gamma_d} N_i V_d\,d\Gamma,
$$

$$
f_i=\int_{\Omega} \rho \Phi_i\,d\Omega .
$$

O uso dos multiplicadores de Lagrange traz algumas desvantagens, pois quebra a propriedade de definida positiva da matriz de rigidez por meio da matriz $G$ e aumenta a ordem do sistema linear a ser resolvido, conforme mostrado no sistema acima. Entretanto, em problemas tridimensionais, o volume cresce a uma taxa muito maior do que a área da fronteira. Assim, pode-se concluir que, à medida que o problema aumenta, a influência dos multiplicadores de Lagrange na ordem do sistema linear tende a se tornar menos significativa.

Como o sistema não é definido positivo, foi utilizado o método GMRES para resolvê-lo [7].

---

## 2.2. Função de forma do EFG

As funções de forma são responsáveis pelas interconexões entre os nós e são fundamentais para a precisão final da solução, pois constituem a base do espaço de solução discretizado. Essas funções são construídas utilizando a técnica dos interpolantes por mínimos quadrados móveis (*moving least squares*) [2].

Ao minimizar uma norma de erro discreta ponderada, obtém-se:

$$
\Phi_I(x)=p^T(x)A^{-1}(x)B_I(x),
$$

em que:

$$
p^T(x)=[1,x,y,z],
$$

$$
A(x)=\sum_{I=1}^{n_x} w(x-x_I)p(x_I)p^T(x_I),
$$

$$
B_I(x)=w(x-x_I)p(x_I).
$$

Como foi escolhido um vetor $p$ de primeira ordem, as funções de forma do EFG também são interpolantes de primeira ordem. Elas possuem suporte compacto devido à natureza da função peso $w(x-x_I)$, que é escolhida de modo a ser diferente de zero apenas em uma pequena parte do domínio.

As funções de forma a serem avaliadas em $x$ são então construídas encontrando-se os nós cujo suporte da função peso contém $x$.

A conectividade dos nós passa a ser considerada porque vários domínios de influência se sobrepõem. A escolha da função peso também é crítica para a precisão da solução. Ela possui alguns parâmetros de escala para controlar seu suporte, que não pode ser tão pequeno a ponto de tornar a matriz $A$ quase singular.

A distribuição dos nós e o tamanho do suporte da função peso, isto é, o tamanho do domínio de influência, devem ser organizados de modo que nenhuma região do domínio permaneça descoberta. Havendo uma distribuição nodal apropriada, esse procedimento desempenha, de maneira simples, um papel semelhante ao refinamento $h$ no FEM.

Neste trabalho, foi utilizada como função peso a função *spline* cúbica:

$$
w(r)=
\begin{cases}
\dfrac{2}{3}-4r^2+4r^3, & 0\le r\le \dfrac{1}{2}, \\
\dfrac{4}{3}-4r+4r^2-\dfrac{4}{3}r^3, & \dfrac{1}{2}<r\le 1, \\
0, & r>1.
\end{cases}
$$

Para fins de comparação, também foi utilizada a função peso quadrática:

$$
w(r)=
\begin{cases}
1-2r^2, & 0\le r\le \dfrac{1}{2}, \\
2(1-r)^2, & \dfrac{1}{2}<r\le 1, \\
0, & r>1.
\end{cases}
$$

Para construir as funções peso em um espaço tridimensional, é necessário escolher se o domínio de influência será esférico ou cúbico. Para construir um domínio de influência esférico, $r$ será a distância euclidiana. Se o domínio cúbico for escolhido, a função de forma final pode ser obtida utilizando o produto tensorial:

$$
w(r_x)w(r_y)w(r_z).
$$

De fato, para definir o suporte da função peso, um parâmetro deve ser colocado no argumento de $w$, como mostrado a seguir. A definição de $d_{m}$ para cada nó ($d_{mI}$) fornece a influência desse nó e, portanto, o suporte de sua função de forma:

$$
r_x=\frac{|x-x_I|}{d_{mI}},
\qquad
r_y=\frac{|y-y_I|}{d_{mI}},
\qquad
r_z=\frac{|z-z_I|}{d_{mI}}.
$$

Outro aspecto importante que deve ser levado em conta é o custo computacional. A construção das funções de forma e o cálculo de suas derivadas constituem a parte mais custosa do EFG em termos de tempo computacional.

Para acelerar o método e encontrar de maneira eficiente os nós no domínio de influência, algumas técnicas de localização de pontos podem ser aplicadas, como mostrado em [8].
