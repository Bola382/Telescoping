# Pastas
## Algoritmo
Arquivos necessários para execução do Telescoping Sampler.
* Funcao Telescoping File.R: Algoritmo TS, utilizar para replicar execuções.
* Relabel File Parallel.R: Soluciona label-switching das replicas.
* Relabeling File.R: Soluciona label-switching.
* Replicando Telescopiong File.R: Controla as replicações do TS.
* Telescoping File.R: Algoritmo TS.
* Verificando Convergencia TS.R: Verificação de convergência através de traceplots e graficos acf.

## Funcoes Auxiliares
Funções necessárias para execução do algoritmo.

## Geracao de dados
Conjunto de dados sintéticos e funções utilizadas para os gerar.

## Outputs
Saídas do algoritmo.

## Outras opcoes para o Telescoping
Iplementações antigas do TS.

# Teoria
## Distribuição t de Student assimétrica
Dizemos que uma variável aleatória $Y$ tem distribuição t de Student assimétrica quando
    $$Y = \mu + X/\sqrt{U}.$$
    Em que $X \sim \mathrm{SN}(0,\sigma^2,\lambda)$ e $U\sim\mathrm{Gama}(\nu/2,\nu/2)$. Escrevemos $Y \sim \mathrm{ST}(\mu,\sigma^2,\lambda,\nu).$
    A v.a. $Y$ tem densidade dada por
    $$\mathrm{St}(y|\mu,\sigma^2,\lambda,\nu) = 2\mathrm{t}(y|\mu,\sigma^2,\nu)\mathrm{T}_{\nu+1}\left(z\lambda \sqrt{\frac{\nu+1}{z^2 + \nu}}\right),$$
    com $z = (y-\mu)/\sigma.$

A distribuição t de Student assimétrica admite a seguinte representação hierárquica:
$$Y|T=t,U=u \sim \mathrm{N}(\mu + \Delta t,u^{-1}\tau^2);$$
$$T|U=u \sim \mathrm{N}_{(0,\infty)}(0,u^{-1});$$
$$U \sim \mathrm{Gama}(\nu/2,\nu/2).$$
Com $\Delta = \sigma\delta$, $\tau^2 = \sigma^2(1-\delta^2)$ e $\delta = \lambda/\sqrt{1+\lambda^2}$.

## Modelo
Seja $Y_i$ a resposta do $i$-ésimo indivíduo, a este indivíduo considere uma variável latente discreta $Z_i$ tal que, dado $Z_i = j$, $Y_i$ dependa do vetor de covariáveis $x_i$ conforme
    $$Y_{i} = x_i^\top\beta_j + \varepsilon_{i},\; j=1,\dots,G.$$
    Onde $G$ é o número de componentes na mistura, $\varepsilon_i|Z_i=j \sim \mathrm{ST}(b\Delta_j,\sigma_j^2,\lambda_j,\nu)$. Supomos graus de liberdade idênticos a todas as componentes por conta da dificuldade de sua estimação. O parâmetro de localização de $\varepsilon_i|Z_i=j$ foi escolhido de forma que o valor esperado dos resíduos fosse zero dada a componente. A constante $b$ é dada por
    $$b = -\sqrt{\frac{\nu}{\pi}}\frac{\Gamma((\nu-1)/2)}{\Gamma(\nu/2)}.$$
As especificações a priori são as seguintes:
    $$G-1 \sim \mathrm{BNB}(1,4,3);$$
    $$\pmb{p}|G \sim \mathrm{Dirichlet}(\gamma/G), \gamma \sim \mathrm{F}(6,3);$$
    $$\beta_j|G \sim \mathrm{N}_{k}(0, 100I_k);$$
    $$\Delta_j|G \sim \mathrm{N}(0,100);$$
    $$\tau_j^2|G \sim \mathrm{InvGama}(0.1,0.1);$$
    $$\nu|\alpha \sim \mathrm{exp}(\alpha), \alpha \sim \mathrm{Unif}(0.02,0.5),$$

O passo 3(a) do algoritmo encontrado em https://doi.org/10.1214/21-BA1294 foi realizado utilizando um algoritmo Gibbs sampler considerando $G$ fixo, construido a partir das prioris mostradas anteriormente juntamente com a verossimilhança do modelo. A condicional completa de $\nu$ necessitou de um passo de MH com proposta lognormal com logmédia no log do passo anterior.

Por ser um modelo de misturas temos problemas de label-switching, podemos permutar os rótulos dos componentes e temos a mesma verossimilhanca, seguimos a sugestão encontrada em https://doi.org/10.1214/21-BA1294.
