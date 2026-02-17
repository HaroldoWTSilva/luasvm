local test = {}
function test.score(trained, test_set)
	local score = 0
	local size = #test_set
	for t=1,size do
		if trained(test_set[t].point) == test_set[t].label then
			score = score + 1
		end
	end
	return score/size
end

test.get_generator_plane =function()
	local ret
	local x = 0.0
	local y = 0.0
	local label
    x = math.random()*10 - 5
    y = math.random()*10 - 5
    label = (x+y > 0.0)
    ret = {point={x,y}, label=label}
	return ret
end

test.get_generator_quadcircle = function()
--O ponto deve estar dentro do quadrado (4,1) (5,8)
--mas fora do circulo de centro (6,4) e raio 1
	local x = 0.0
	local y = 0.0
	local rx
	local ry
	local point
    x = math.random()*20 - 10
    y = math.random()*20 - 10
    rx = x - 6
    ry = y - 3
    point = {point={x/10.0, y/10.0}}
    if x > 4.0 and x < 8.0 and y > 1.0 and y < 5
        and (rx * rx + ry * ry) > 1.0 then
        point.label = true
    else
        point.label = false
    end
    return point
end


test.run = function(point_generator, training_size, fitfunc)
    local training_set = {}
    local test_set = {}

    for i=1,training_size do
        training_set[i] = point_generator()
        test_set[i] = point_generator()
    end
    print("Treinando com", training_size, "pontos!")

    local trained = fitfunc(training_set)

    print("Testando no conjunto de treino...")
    print("Pontuação: ", test.score(trained, training_set) )
    print("Testando no conjunto de teste...")
    print("Pontuação: ", test.score(trained, test_set) )
end

return test
