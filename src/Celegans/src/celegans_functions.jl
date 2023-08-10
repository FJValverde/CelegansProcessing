# This part of the module gathers all transformation functions, etc.
# FVA: Probably this should be done 
    

"""
        connectomeToMultivaluedMatrix(dfConnectome::DataFrame)

A function to create two matrices for storing the connectivity
information of a connectome in a DataFrame `dfConnectome`.

See also [`connectomeToBinaryMatrix`](@ref)
"""
function connectomeToMultivaluedMatrix(dfConnectome::DataFrame)
    mat_num = zeros(nNeurons, nNeurons)
    mat_bin = zeros(nNeurons, nNeurons)
    for link in eachrow(dfConnectome)
        from_index = link[:IndexSending]
        to_index = link[:IndexReceiving]
        mat_num[from_index, to_index] = link[:Number]
        mat_bin[from_index, to_index] = 1.0        
    end
    return(mat_num,mat_bin)
end


"""
        connectomeToBinaryMatrix(dfConnectome::DataFrame)

A function to create a binary matrix for storing the connectivity
information of a connectome in a DataFrame `dfConnectome`.

See also [`connectomeToMultivaluedMatrix`](@ref)
"""
function connectomeToBinaryMatrix(dfConnectome::DataFrame)
    mat_bin = zeros(nNeurons, nNeurons)
    for link in eachrow(dfConnectome)
        from_index = link[:IndexSending]
        to_index = link[:IndexReceiving]
         mat_bin[from_index, to_index] = 1.0        
    end
    return(mat_bin)
end
