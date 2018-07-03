useGPU = true;
dataType = 'single';
backward = true;
rec1 = RecurrentLayer(struct('hidden_dim',512,'input_dim',512,'useGPU',useGPU,'dataType',dataType,'backward',backward));
rec2 = RecurrentLayer(struct('hidden_dim',512,'input_dim',512,'useGPU',useGPU,'dataType',dataType,'backward',backward));
