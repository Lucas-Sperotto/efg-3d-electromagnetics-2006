# 4. Conclusão

Neste trabalho, o método sem malha **Element-Free Galerkin** (**EFG**) foi aplicado para obter a solução numérica de problemas eletrostáticos tridimensionais (3-D).

Foi apresentada a formulação com multiplicadores de Lagrange para impor a condição de contorno de Dirichlet. Para ilustrar a efetividade do método, foram apresentados resultados numéricos para um problema eletrostático.

A norma do erro dos valores do potencial, para um problema típico, foi de **1,40%** com o EFG, enquanto o FEM apresentou erro de **3,56%** utilizando os mesmos nós.

A formulação EFG envolve muitos parâmetros, como a função peso, o domínio de influência e a integração numérica. O controle desses parâmetros é uma questão fundamental para a implementação efetiva do método. Neste artigo, foram estudados os efeitos desses parâmetros em domínios tridimensionais.

Foi aplicada uma estratégia para utilizar diferentes tamanhos de domínio de influência de acordo com a densidade local dos nós, assegurando uma construção eficiente e correta das funções de forma.

A integração da forma fraca também foi estudada, e foram apresentados resultados mostrando como ocorre sua convergência.

Esses resultados ampliam o conhecimento sobre como e onde os parâmetros do EFG influenciam a precisão da solução e facilitam a aplicação do método EFG em domínios tridimensionais com nuvens gerais de nós.

---

## Agradecimentos

Este trabalho recebeu apoio parcial das agências brasileiras **CAPES** e **CNPq**.

---

## Referências

[1] T. Belytschko, Y. Krongauz, D. Organ, M. Fleming e P. Krysl, “Meshless methods: An overview and recent developments,” *Computer Methods in Applied Mechanics and Engineering*, v. 139, n. 1–4, p. 3–47, maio 1996.

[2] T. Belytschko, Y. Y. Lu e L. Gu, “Element-free Galerkin methods,” *International Journal for Numerical Methods in Engineering*, v. 37, p. 229–256, 1994.

[3] S. Suleau e P. Bouillard, “One-dimensional dispersion analysis for the element-free Galerkin method for the Helmholtz equation,” *International Journal for Numerical Methods in Engineering*, v. 47, n. 6, p. 1169–1188, fevereiro 2000.

[4] T. Belytschko, Y. Lu e L. Gu, “Crack propagation by element-free Galerkin methods,” *Engineering Fracture Mechanics*, v. 51, n. 2, p. 295–315, maio 1995.

[5] V. Cingoski, N. Miyamoto e H. Yamashita, “Element-free Galerkin method for electromagnetic field computations,” *IEEE Transactions on Magnetics*, v. 34, n. 5, p. 3236–3239, setembro 1998.

[6] S. Fernández-Méndez e A. Huerta, “Imposing essential boundary conditions in mesh-free methods,” *Computer Methods in Applied Mechanics and Engineering*, v. 193, n. 12–14, p. 1257–1275, março 2004.

[7] Y. Saad e M. Schultz, “GMRES: A generalized minimal residual algorithm for solving nonsymmetric linear systems,” *SIAM Journal on Scientific Computing*, v. 7, n. 3, p. 856–869, julho 1986.

[8] G. F. Parreira, A. R. Fonseca, A. C. Lisboa, E. J. Silva e R. C. Mesquita, “Efficient algorithms and data structures for element-free Galerkin method,” em *Proceedings of the 15th Conference on the Computation of Electromagnetic Fields*, junho 2005, p. 212–213.
