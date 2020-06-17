feria: feria.h feria_parameter.o feria_mesh.o feria_env.o feria_sky.o feria_PV.o feria_fitsio.o feria_main.cpp
	g++ -Wall -O2 -L/opt/local/lib/ -I/opt/local/include/ -lcfitsio feria_parameter.o feria_mesh.o feria_env.o feria_sky.o feria_PV.o feria_fitsio.o feria_main.cpp -o feria

feria_fitsio.o: feria.h feria_parameter.o feria_sky.o feria_fitsio.cpp
	g++ -Wall -O2 -I/opt/local/include/ -c feria_fitsio.cpp

feria_PV.o: feria.h feria_parameter.o feria_PV.cpp
	g++ -Wall -O2 -I/opt/local/include/ -c feria_PV.cpp

feria_sky.o: feria.h feria_parameter.o feria_sky.cpp
	g++ -Wall -O2 -I/opt/local/include/ -c feria_sky.cpp

feria_env.o: feria.h feria_parameter.o feria_env.cpp
	g++ -Wall -O2 -I/opt/local/include/ -c feria_env.cpp

feria_mesh.o: feria.h feria_parameter.o feria_mesh.cpp
	g++ -Wall -O2 -I/opt/local/include/ -c feria_mesh.cpp

feria_parameter.o: feria.h feria_parameter.cpp
	g++ -Wall -O2 -I/opt/local/include/ -c feria_parameter.cpp

clean:
	rm -f *.o feria

