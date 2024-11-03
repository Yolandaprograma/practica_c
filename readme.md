## OPERACIONS AMB MATRIUS I VECTORS EN C

#### L'objectiu de la pràctica serà realitzar una biblioteca que permeti realitzar operacions matricials/vectorials utilitzant el llenguatge de programació C. La llibreria haurà de permetre realitzar les operacions següents:
- Multiplicació escalar: multiplicar un vector per un escalar.
- Producte escalar: producte punt entre dos vectors.
- Magnitud: calcular la magnitud d’un vector.
- Ortogonals: determinar si dos vectors són ortogonals.
- Projecció: calcular la projecció d’un vector sobre un altre.
- Norma de Frobenius: calcular la norma de Frobenius d’una matriu.
- Infini-norma: calcular la infini-norma d’una matriu.
- Norma-u: calcular la norma ú d’una matriu.
- Diagonal dominant: determinar si una matriu és o no diagonal dominant.
- Multiplicació matriu per vector.
- Resoldre sistemes d’equacions lineals: amb un mètode iteratiu.

#### Com estem fent el desenvolupament d’aquesta biblioteca, cal tenir un conjunt de proves per validar la funcionalitat de cada funció. Per fer això, a part de les funcions anteriors, crearem un conjunt de matrius i vectors i una funció d’inicialització per aquestes estructures de dades.
#### Les estructures i constants que declarareu són les següents:
- Constant N igual a 512.
- Matrius (Mat i MatDD) de números reals (punt flotant de precisió simple) de NxN elements.
- Vectors (V1, V2, V3 i V4) de números reals de N elements.
#### Amb l'objectiu que tots els alumnes utilitzeu les mateixes dades d'entrada, a l'hora de generar els elements que contindran la matriu i els vectors, la funció initData haurà d'utilitzar les funcions incloses a l'estàndard de C rand() i srand() per generar els valors dels elements. La funció rand() permet generar números aleatoris utilitzant un generador de nombres aleatoris. La funció srand() permet instanciar la llavor que utilitzarà com a base el generador de nombre aleatoris. Per això, passarem com a entrada la llavor “334411” (srand(334411)). Amb això aconseguirem generar números pseudo-aleatoris, que seguiran la mateixa seqüència en cada execució del programa. Aquesta funció requereix ser cridada abans de la funció rand().
#### La funció d'inicialització la proporcionem ja programada (tots heu de fer servir aquesta funció per poder comprovar que els resultats són correctes):
- void PrintVect( float vect[N], int from, int numel )
- void PrintRow( float mat[N][N], int row, int from, int numel )
- void MultEscalar( float vect[N], float vectres[N], float alfa )
- float Scalar( float vect1[N], float vect2[N] )
- float Magnitude( float vect[N] )
- int Ortogonal( float vect1[N], float vect2[N] )
- void Projection( float vect1[N], float vect2[N], float vectres[N] )
- float Infininorm( float M[N][N] )
- float Onenorm( float M[N][N] )
- float NormFrobenius( float M[N][N] )
- int DiagonalDom( float M[N][N] )
- void Matriu_x_Vector( float M[N][N], float vect[N], float vectres[N] )
- int Jacobi( float M[N][N] , float vect[N], float vectres[N], unsigned iter )
