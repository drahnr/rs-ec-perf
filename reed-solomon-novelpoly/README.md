# reed-solomon-novelpoly

An implementation of  Novel Polynomial Basis and its Application to Reed-Solomon Erasure Codes [1] [2] .

Runs encoding and reconstruction in `O(n lg(n))`. Note that for small number `n` there is a static offset due to a walsh transform over the full domain in reconstruction.

## Goals

Be fast for `n > 100` and in particular for `100 < n < 2_000`.

[1]: https://www.citi.sinica.edu.tw/papers/whc/4454-F.pdf
[2]: https://arxiv.org/abs/1404.3458

# License

Apache-2.0 and MIT
