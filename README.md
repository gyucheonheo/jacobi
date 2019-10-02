# JACOBI BENCHMARKING

## HOW TO COMPILE

```bash
$>make jacobi    /* Serial Version */
$>make jacobi_mp /* openMP Version */
```

## HOW TO RUN

```bash
$>./jacobi -n <size> -i <iterations> -c <convergence> /* Serial Version */
$>./jacobi_mp -t <threads> -n <size> -i <iterations> -c <convergence> /* OpenMP Version */
```

## HOW TO TEST

```bash
$>make test    /* Serial Version */
$>make test_mp /* OpenMP Version */
```



