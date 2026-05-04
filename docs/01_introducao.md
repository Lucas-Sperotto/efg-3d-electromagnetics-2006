**Navegação:** [README](../README.md) | [Capa e resumo](00_capa_resumo.md) | [1. Introdução](01_introducao.md) | [2. Formulação matemática](02_formulacao_matematica.md) | [3. Resultados numéricos](03_resultados_numericos.md) | [4. Conclusão e referências](04_conclusao_referencias.md)

---

# 1. Introdução

Problemas eletromagnéticos em domínios tridimensionais (3-D) são muito mais complexos do que em dimensões inferiores, devido ao aumento do custo computacional e às dificuldades na construção da discretização do domínio.

O método dos elementos finitos (**FEM**, do inglês *Finite Element Method*) apresenta características que permitem considerar diferentes tipos de materiais e modelar domínios geometricamente complexos. Entretanto, uma malha de boa qualidade é necessária para produzir resultados precisos. Esse aspecto se torna ainda mais crítico na solução de problemas com múltiplas escalas de comprimento, como em ensaios não destrutivos por micro-ondas, ou em problemas nos quais o remalhamento é necessário, como aqueles que envolvem fronteiras móveis.

Nesse contexto, os métodos sem malha ampliam a área de aplicação do eletromagnetismo computacional. Esses métodos fornecem uma solução utilizando apenas informações nodais, sem a necessidade de uma malha para estabelecer a conectividade entre os nós.

Existem diversos métodos sem malha, tais como o método das partículas com núcleo reprodutor (**RKPM**, *Reproducing Kernel Particle Method*), o método dos elementos difusos (**DEM**, *Diffuse-Element Method*) e o método da partição da unidade (**PUM**, *Partition of Unity Method*), cada um com características próprias. Uma revisão desses métodos pode ser encontrada em [1].

Entre esses métodos, um dos mais promissores é o método de Galerkin sem elementos (**EFG**, *Element-Free Galerkin*) [2], devido às boas taxas de convergência, à facilidade de criar a discretização e representar o domínio, além da independência da integração da forma fraca em relação a uma malha de conectividade.

Em trabalhos anteriores, o EFG foi aplicado com sucesso nas áreas de acústica [3], mecânica [4] e eletromagnetismo [5]. Entretanto, esses trabalhos apresentam implementações apenas para problemas bidimensionais e não indicam se seus resultados podem ser estendidos para três dimensões.

O objetivo deste artigo é aplicar o método EFG à solução de problemas eletromagnéticos em domínios tridimensionais (3-D). Em particular, o EFG é aplicado à solução de um problema eletrostático envolvendo múltiplos cantos. Esse problema é utilizado como referência (*benchmark*) para definir diversos parâmetros do EFG, tais como a função peso, o tamanho do domínio de influência e o número de células de integração, os quais não são fáceis de escolher e exercem forte influência sobre a precisão da solução.

Para atingir esse objetivo, primeiro é apresentada a formulação matemática, introduzindo as formas fraca e discreta do problema na Seção II. As funções de forma do EFG também são discutidas nessa seção. A Seção III apresenta uma técnica para definir o suporte das funções de forma com base na densidade local dos nós e discute resultados numéricos relacionados aos parâmetros do EFG.

---

**Anterior:** [Capa e resumo](00_capa_resumo.md) | **Próximo:** [2. Formulação matemática](02_formulacao_matematica.md)
