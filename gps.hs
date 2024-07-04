{-Funcion main-}

main :: IO ()
main = do
  texto <- readFile "sat.txt"
  let cartesianas = kmAm (iterar5veces 0 [0, 0, 0] (iniciarSistema texto)) -- Se obtienen las coordenadas cartesianas
      geodesicas = cartAgeo (primeraIter cartesianas) cartesianas 0 -- Se obtienen las coordenadas geodesicas
      geodesicasGrad = map radianesAgrados geodesicas
      s0 = "\n---------Posicion de un GPS a partir de las coordenadas y el tiempo de 10 satélites GPS ---------"
      s = s0++"\nCartesianas (X,Y,Z) en metros = " ++ show cartesianas ++ "\nGeodesicas (Long, Lat, H) en grados, H en metros = " ++ show [(geodesicasGrad!!0), (geodesicasGrad!!1), (geodesicas!!2)]
      s2 = s++"\nUTM (E,N,Zona) = " ++ show (utm (radianesAgrados(geodesicas!!0)) (geodesicas!!0) (geodesicas!!1) )
  putStrLn s2

gradosAradianes :: Double -> Double
gradosAradianes val = val * (pi/ 180) 

radianesAgrados :: Double -> Double
radianesAgrados val = val * (180/pi)

mAkm :: Double -> Double
mAkm x = x*1000

iniciarSistema :: String -> [[Double]]
iniciarSistema texto = sis
  where mat = leerSat texto
        datosIni = iniSist mat 
        sis = moverRoi datosIni --Es necesario pasar los roi a la derecha para el funcionamiento

moverRoi :: [[Double]] -> [[Double]]
moverRoi sis = map rotaFila sis 
rotaFila :: [Double] -> [Double]
rotaFila [p,x,y,z] = [x,y,z,p]

{- Funciones para leer un fichero y convertirlo en lista de doubles-}
leerSat :: String -> [[String]]
leerSat texto = matriz 
  where lineas = lines texto
        matriz = map words lineas

iniSist :: [[String]] -> [[Double]]
iniSist strSist = map listaAfloat strSist

listaAfloat :: [String] -> [Double]
listaAfloat str = map toDouble str

toDouble :: String -> Double
toDouble str = read str
{- Convierte un entero a Coma Flotante-}
intAdouble :: Int -> Double
intAdouble x = fromIntegral x

{-------------------------------------------------------------------}

{- Funciones para convertir de WGS84 a UTM -}

utm :: Double -> Double -> Double -> (Double, Double, Int)
utm longg long lat = (este,norte,zona)
  where t = sinh( atanh(sin(lat)) - (((2 * sqrt(n)) / (1+n)) * atanh(((2 * sqrt(n)) / (1+n))*sin(lat) )))
        xi' = atan(t/(cos(long-long0)))
        eta' = atanh (sin(long-long0) / sqrt(1 + (t^2)))
        sigma = 1 + (sigma1 xi' eta' alfa 0 1)
        tau = tau1 xi' eta' alfa 0 1
        e = (2 * sqrt(n)) / (1+n)
        long0 = gradosAradianes long0gr
        long0gr = (intAdouble zona) * 6.0 - 183.0 {-se calcula la longitud inicial-}
        zona = floor ((longg+186)/6) 
        este = 500000 + k0*aMayus*(eta'+(este1 xi' eta' alfa 0.0 1))   --en metros
        norte = (determinaN0 lat) + k0*aMayus*(xi'+(norte1 xi' eta' alfa 0.0 1))
        alfa = [alfa1,alfa2,alfa3]
        alfa1 = (n/2) - ((2*(n^2))/3) + ((5*(n^3))/16)
        alfa2 = ((13*(n^2))/48) - ((3*(n^3))/5)
        alfa3 = (61*(n^3))/240
        f = 0.0033528106647474
        n = (f / (2 - f))
        k0 = 0.9996
        aMayus = (aMinus/(1+n))*(1+((n^2)/4)+((n^4)/64))
        aMinus = 6378137 --en metros

este1 :: Double -> Double -> [Double] -> Double -> Int -> Double
este1 xi' eta' alfa acc j = if j==4 then acc
                           else este1 xi' eta' alfa acc1 (j+1)
  where acc1 = acc + (alfa!!(j-1))*cos(2*j2*xi')*sinh(2*j2*eta')
        j2 = intAdouble j

norte1 :: Double -> Double -> [Double] -> Double -> Int -> Double
norte1 xi' eta' alfa acc j = if j==4 then acc
                           else norte1 xi' eta' alfa acc1 (j+1)
  where acc1 = acc + (alfa!!(j-1))*sin(2*j2*xi')*cosh(2*j2*eta') 
        j2 = intAdouble j

sigma1 :: Double -> Double -> [Double] -> Double -> Int -> Double
sigma1 xi' eta' alfa acc j = if j==4 then acc
                           else sigma1 xi' eta' alfa acc1 (j+1)
  where acc1 = acc + (2*j2*(alfa!!(j-1))*cos(2*j2*xi')*cos(2*j2*eta'))
        j2 = intAdouble j

tau1 :: Double -> Double -> [Double] -> Double -> Int -> Double
tau1 xi' eta' alfa acc j = if j==4 then acc
                         else tau1 xi' eta' alfa acc1 (j+1)
  where acc1 = acc + (2*j2*(alfa!!(j-1))*sin(2*j2*xi')*sinh(2*j2*eta'))
        j2 = intAdouble j

determinaN0 :: Double -> Double
determinaN0 lat = if lat>=0 then 0 else 10000000 --en metros

{--------------------------------------------------------------------}

{- Funciones para convertir de XYZ a WGS84 -}

kmAm :: [Double] -> [Double] 
kmAm [x,y,z] = [(x*1000),(y*1000),(z*1000)] 

primeraIter :: [Double] -> [Double]
primeraIter [x,y,z] = [long, lat, h, psi]
  where p = sqrt((x^2) + (y^2))
        lat = atan (z/p)
        long = atan2 y x
        h=0
        re=6378137
        rp=6356752.314
        psi=atan((re*(z-(h*sin(lat))))/(rp*(p-(h*cos(lat)))))

cartAgeo :: [Double] -> [Double] -> Int -> [Double]
cartAgeo [long, lat, h, psi] [x,y,z] i = if i==4 then [long, lat, h]
                                         else cartAgeo [long, lat1, h1, psi1] [x,y,z] (i+1)
  where re=6378137
        rp=6356752.314
        p = sqrt((x^2) + (y^2))
        lat1= atan((re*tan(psi))/rp)
        h1 = (p-(re*cos(psi)))/(cos(lat1))
        psi1 = atan(re*((z-(h1*sin(lat1)))/(rp*(p-(h1*cos(lat1))))))
{--------------------------------------------------------------------}


{- Funciones para resolver los sumatorios y formar un sistema de ecuaciones -}
distancia :: [Double] -> [Double] -> Double
distancia coord0 coordi = sqrt(d1 + d2 + d3)
  where d1 = ((coord0!!0) - (coordi!!0))^2
        d2 = ((coord0!!1) - (coordi!!1))^2
        d3 = ((coord0!!2) - (coordi!!2))^2

-- Itera 5 veces sobre creaSistema aproximando (x,y,z)
iterar5veces :: Int -> [Double] -> [[Double]] -> [Double] 
iterar5veces i coord0 coordi = if i==5 then coord0 else iterar5veces (i+1) coordsig coordi
                               where sistema0 = [[0,0,0,0],[0,0,0,0],[0,0,0,0]]
                                     coordsig = sumarFilas coord0 (extraeSolucion (gauss (creaSistema1 coord0 coordi sistema0 0) 0))

-- Extrae el valor de (x,y,z) de la matriz reducida
extraeSolucion :: [[Double]] -> [Double]
extraeSolucion sistema = [((sistema!!0)!!3),((sistema!!1)!!3),((sistema!!2)!!3)]

-- Calcula y acumula el sumatorio de los coeficientes del sistema de ecuaciones para los 10 satelites
creaSistema1 :: [Double] -> [[Double]] -> [[Double]] -> Int -> [[Double]] 
creaSistema1 coord0 coordi sistema i | i==10 = sistema
                                     | otherwise = creaSistema1 coord0 coordi acc (i+1) 
                                      where acc = sumarSistemas sistema (creaSistema2 coord0 (coordi!!i) ((coordi!!i)!!3))

-- Calcula los coeficientes de un satelite para el proximo sumatorio
creaSistema2 :: [Double] -> [Double] -> Double -> [[Double]]
creaSistema2 [x0,y1,z1] [x,y,z] roi = [[x1,y1,z1,r1],[x2,y2,z2,r2],[x3,y3,z3,r3]]
  where d = distancia coord0 coordi
        x1 = (((x0)-(x))^2)/(d^2)
        y1 = (((y0)-(y))*((x0)-(x)))/(d^2)
        z1 = (((z0)-(y))*((x0)-(x)))/(d^2)
        r1 = ((roi - d)*((x0)-(x)))/d
        x2 = (((x0)-(x))*((y0)-(y)))/(d^2)
        y2 = (((y0)-(y))^2)/(d^2)
        z2 = (((z0)-(y))*((y0)-(y)))/(d^2)
        r2 = ((roi - d)*((y0)-(y)))/d
        x3 = (((x0)-(x))*((z0)-(y)))/(d^2)
        y3 = (((y0)-(y))*((z0)-(y)))/(d^2)
        z3 = (((z0)-(y))^2)/(d^2)
        r3 = ((roi - d)*((z0)-(y)))/d

{-Funcion que suma dos sistemas-}
sumarSistemas :: [[Double]] -> [[Double]] -> [[Double]] 
sumarSistemas s1 s2 = [suma1, suma2, suma3]
  where suma1 = sumarFilas (s1!!0) (s2!!0)
        suma2 = sumarFilas (s1!!1) (s2!!1)
        suma3 = sumarFilas (s1!!2) (s2!!2)

{- Funciones para resolver el sistema de ecuaciones por Gauss -}
--Dividir todos los elementos de una lista entre un valor
dividirFila :: [Double] -> Double -> [Double]
dividirFila [] _  = []
dividirFila (x:xs) y = (x/y):(dividirFila xs y)
--Multiplica una fila por un Double
multiplicarFila :: [Double] -> Double -> [Double] 
multiplicarFila [] _  = []
multiplicarFila (x:xs) y = (x*y):(multiplicarFila xs y)
--Restar dos filas
restarFilas [] [] = []
restarFilas (x:xs) (y:ys) = (x - y):(restarFilas xs ys)
--Sumar dos filas
sumarFilas [] [] = []
sumarFilas (x:xs) (y:ys) = (x + y):(sumarFilas xs ys)

-- funcion que dada una fila que se desea reducir y una fila con pivote 1 opera para reducir la fila deseada 
reducirFila :: [Double] -> [Double] -> Double -> [Double]
reducirFila filaAreducir filaPivote coeficiente = restarFilas filaAreducir filaMultiplicada
  where filaMultiplicada = multiplicarFila filaPivote coeficiente

{-Funcion que recibe una lista de tres listas y las rota 1 posicion-}
intercambiaFilas :: [[Double]] -> [[Double]] 
intercambiaFilas sistema = [(sistema!!2),(sistema!!0),(sistema!!1)]

{-Funcion que recibe el sistema y la columna que se desea reducir y gira el sistema en torno a un pivote
que coincida con la fila-}
recoloca :: [[Double]] -> Int -> [[Double]]
recoloca sistema columna = if ((sistema!!columna)!!columna)==0 then recoloca (intercambiaFilas sistema) columna 
                           else sistema 

--Funcion que resuelve un sistema de tres ecuaciones y tres incognitas por gauss-jordan
gauss :: [[Double]] -> Int -> [[Double]]
gauss sistema 3 = sistema
gauss sistema columna | columna==0 = gauss [d, r1, r2] c1
                      | columna==1 = gauss [r2, d, r1] c1
                      | otherwise = gauss [r1, r2, d] c1
                      where sis = recoloca sistema columna
                            d = dividirFila (sis!!columna) ((sis!!columna)!!columna)  
                            r1 = reducirFila (sis!!((columna+1)`mod`3)) d ((sis!!((columna+1)`mod`3))!!columna) 
                            r2 = reducirFila (sis!!((columna+2)`mod`3)) d ((sis!!((columna+2)`mod`3))!!columna)
                            c1 = columna + 1

{--------------------------------------------------------------------}