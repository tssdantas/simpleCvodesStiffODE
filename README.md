my simple CVODES example
===============

Um exemplo simples que compila e resolve um problema de uma ODE stiff com o CVODES


## Instalação do CVODES
No terminal do WSL ou Linux (Ubuntu 22.04) inserir os seguintos comandos
* sudo apt-get update
* sudo apt-get install libsundials-dev 
* (se não funcionar tentar sudo apt-get install libsundials-cvodes4)


## Instruções para compilar o exemplo
No terminal do WSL ou Linux inserir os seguintos comandos
* make
* ./stiff_ode
