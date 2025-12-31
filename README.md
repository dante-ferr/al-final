# Extrator de Kernel de Matrizes

Ferramenta em C para calcular a dimensão e a base do núcleo (kernel) de uma matriz $m \times n$, utilizando eliminação de Gauss-Jordan para obter a forma escalonada reduzida por linhas (RREF).

Desenvolvido para a disciplina de Álgebra Linear.

## Requisitos

- **Compilador C** (GCC ou Clang)
- **CMake** (Versão 3.10 ou superior)

## Compilação

Este projeto utiliza o CMake para gerar os arquivos de build. Na raiz do projeto, execute:

```bash
cmake -B build
cmake --build build
```

## Execução

O binário será gerado dentro da pasta `build`.

```bash
./build/main
```

### Formato de Entrada

O programa espera a entrada no seguinte formato via `stdin` (pode ser digitado ou colado):

1. **Dimensões:** Dois inteiros representando linhas e colunas.
2. **Matriz:** Os elementos da matriz (floats/doubles), linha por linha.

#### Exemplo de Uso

Para uma matriz $3 \times 3$:

```text
3 3
1 2 3
4 5 6
7 8 9
```

**Saída:**

```text
[Entrada]
|    1.000    2.000    3.000 |
|    4.000    5.000    6.000 |
|    7.000    8.000    9.000 |

[RREF]
|    1.000    0.000   -1.000 |
|    0.000    1.000    2.000 |
|    0.000    0.000    0.000 |

Dimensao Kernel: 1
[Base]
|    1.000   -2.000    1.000 |
```
