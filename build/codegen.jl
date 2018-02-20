using IMEXRKCB

# Uncomment appropriately below to see the code generated. When strip=false
# the code is not simplified and all operations of the implementations are
# generated

embed = false
# impl = IMEXRK3R2R(IMEXRKCB2, embed)
impl = IMEXRK3R2R(IMEXRKCB3e, embed)
# impl = IMEXRK3R2R(IMEXRKCB3c, embed)
# impl = IMEXRK4R3R(IMEXRKCB4,  embed)

@show codegen(impl; strip=true)