**Navegação:** [README](../README.md) | [Capa e resumo](00_capa_resumo.md) | [1. Introdução](01_introducao.md) | [2. Formulação matemática](02_formulacao_matematica.md) | [3. Resultados numéricos](03_resultados_numericos.md) | [4. Conclusão e referências](04_conclusao_referencias.md)

---

# 4. Conclusão

Neste trabalho, o método sem malha **Element-Free Galerkin** (**EFG**) foi aplicado para obter a solução numérica de problemas eletrostáticos tridimensionais (3-D).

Foi apresentada a formulação com multiplicadores de Lagrange para impor a condição de contorno de Dirichlet. Para ilustrar a efetividade do método, foram apresentados resultados numéricos para um problema eletrostático.

No problema de referência, a norma relativa do erro dos valores do potencial foi de **1,40%** com o EFG, enquanto o FEM apresentou erro de **3,56%** utilizando os mesmos nós.

A formulação EFG envolve muitos parâmetros, como a função peso, o domínio de influência e a integração numérica. O controle desses parâmetros é uma questão fundamental para a implementação efetiva do método. Neste artigo, foram estudados os efeitos desses parâmetros em domínios tridimensionais.

Foi aplicada uma estratégia para utilizar diferentes tamanhos do domínio de influência de acordo com a densidade local dos nós, assegurando uma construção eficiente e correta das funções de forma.

A integração da forma fraca também foi estudada, e foram apresentados resultados mostrando como ocorre sua convergência.

Esses resultados ampliam o conhecimento sobre como e onde os parâmetros do EFG influenciam a precisão da solução e facilitam a aplicação do método EFG em domínios tridimensionais com nuvens gerais de nós.

---

## Agradecimentos

Este trabalho recebeu apoio parcial das agências brasileiras **CAPES** e **CNPq**.

---

## Referências

[1] T. Belytschko, Y. Krongauz, D. Organ, M. Fleming e P. Krysl, “Meshless methods: An overview and recent developments,” *Computer Methods in Applied Mechanics and Engineering*, v. 139, n. 1–4, p. 3–47, maio 1996. DOI: [10.1016/S0045-7825(96)01078-X](<https://doi.org/10.1016/S0045-7825(96)01078-X>).

[2] T. Belytschko, Y. Y. Lu e L. Gu, “Element-free Galerkin methods,” *International Journal for Numerical Methods in Engineering*, v. 37, p. 229–256, 1994. DOI: [10.1002/nme.1620370205](https://doi.org/10.1002/nme.1620370205).

[3] S. Suleau e P. Bouillard, “One-dimensional dispersion analysis for the element-free Galerkin method for the Helmholtz equation,” *International Journal for Numerical Methods in Engineering*, v. 47, n. 6, p. 1169–1188, fevereiro 2000. DOI: [10.1002/(SICI)1097-0207(20000228)47:6&lt;1169::AID-NME824&gt;3.0.CO;2-9](<https://doi.org/10.1002/(SICI)1097-0207(20000228)47:6%3C1169::AID-NME824%3E3.0.CO;2-9>).

[4] T. Belytschko, Y. Lu e L. Gu, “Crack propagation by element-free Galerkin methods,” *Engineering Fracture Mechanics*, v. 51, n. 2, p. 295–315, maio 1995. DOI: [10.1016/0013-7944(94)00153-9](<https://doi.org/10.1016/0013-7944(94)00153-9>).

[5] V. Cingoski, N. Miyamoto e H. Yamashita, “Element-free Galerkin method for electromagnetic field computations,” *IEEE Transactions on Magnetics*, v. 34, n. 5, p. 3236–3239, setembro 1998. DOI: [10.1109/20.717759](https://doi.org/10.1109/20.717759).

[6] S. Fernández-Méndez e A. Huerta, “Imposing essential boundary conditions in mesh-free methods,” *Computer Methods in Applied Mechanics and Engineering*, v. 193, n. 12–14, p. 1257–1275, março 2004. DOI: [10.1016/j.cma.2003.12.019](https://doi.org/10.1016/j.cma.2003.12.019).

[7] Y. Saad e M. Schultz, “GMRES: A generalized minimal residual algorithm for solving nonsymmetric linear systems,” *SIAM Journal on Scientific and Statistical Computing*, v. 7, n. 3, p. 856–869, julho 1986. DOI: [10.1137/0907058](https://doi.org/10.1137/0907058).

[8] G. F. Parreira, A. R. Fonseca, A. C. Lisboa, E. J. Silva e R. C. Mesquita, “Efficient algorithms and data structures for element-free Galerkin method,” em *Proceedings of the 15th Conference on the Computation of Electromagnetic Fields*, junho 2005, p. 212–213. DOI: não localizado para a versão em anais.

[8a] G. F. Parreira, A. R. Fonseca, A. C. Lisboa, E. J. Silva e R. C. Mesquita, “Efficient algorithms and data structures for element-free Galerkin method,” *IEEE Transactions on Magnetics*, v. 42, n. 4, p. 659–662, abril 2006. DOI: [10.1109/TMAG.2006.871432](https://doi.org/10.1109/TMAG.2006.871432).

---

**Anterior:** [3. Resultados numéricos](03_resultados_numericos.md) | **Voltar ao:** [README](../README.md)
