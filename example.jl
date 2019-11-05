example_qaplib_format = """
3

1 2 5
2 0 1
4 1 2

1 0 1
0 1 0
1 0 1
"""

# Load problem utilities
include("./src/problems/QAP/qap.jl")

# Parse instance
instance = parse_qap_qaplib(Int64, example_qaplib_format)

# Load approach(es)
include("./src/approaches/pGOMEA/pGOMEA.jl")
include("./src/approaches/qGOMEA/qGOMEA.jl")
include("./src/approaches/SimpleGA/IPSimpleGA.jl")
include("./src/approaches/SimpleGA/RKSimpleGA.jl")

# Find solutions!
optimize_pgomea(bb_wrap_qap(instance), instance.n, 1.0)
optimize_qgomea(bb_wrap_qap(instance), instance.n, 1.0)
optimize_rksimplega(bb_wrap_qap(instance), instance.n, 1.0)
optimize_ipsimplega(PMX(instance.n), bb_wrap_qap(instance), instance.n, 1.0)
optimize_ipsimplega(CX(instance.n), bb_wrap_qap(instance), instance.n, 1.0)
optimize_ipsimplega(OX(instance.n), bb_wrap_qap(instance), instance.n, 1.0)
optimize_ipsimplega(ER(instance.n), bb_wrap_qap(instance), instance.n, 1.0)
