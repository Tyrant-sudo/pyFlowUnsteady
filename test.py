from julia import Julia
jl = Julia(compiled_modules=False)

# 导入 Julia 模块
vlm = jl.eval("import FLOWVLM as vlm; vlm")
read_polar = vlm.ap.read_polar
