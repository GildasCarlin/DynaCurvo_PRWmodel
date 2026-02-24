/**
 * File              : maths_utils.js
 * Original article  : doi
 * Authors           : Gildas Carlin, Ian Manifacier and Jean-Louis Milan
 * Date              : 11/12/2025
 * Description       : Ce fichier regroupe plusieurs fonctions mathématiques utilisées dans le code
 */
//

function getRandomVector(){   
  // Box-Muller transform 
  let u = 0.0, v = 0.0; 
  while(u === 0.0) u = random();  // Converting [0,1) to (0,1)
  while(v === 0.0) v = random();
  let noise = createVector(
    sqrt(-2.0*log(u)) * cos(2.0*PI*v),
    sqrt(-2.0*log(u)) * sin(2.0*PI*v)
  );
  return noise;
}


function NormalizeVector(vecs){
    // takes as input: a vector or an array of vectors
    // returns: the same vector(s) normalized (unit-length)
    // note: if the norm of a vector is almost zero, the null vector is returned

  const isArray=Array.isArray(vecs);
  const array_of_vectors=isArray ? vecs : [vecs];

  let normalization_done=[];
  for (let i=0; i<array_of_vectors.length; i++){
    if (array_of_vectors[i].mag() < pow(10.0, -10)){
      normalization_done.push(createVector(0.0, 0.0));
    } else {
      normalization_done.push(array_of_vectors[i].normalize());
    }
  }
  return isArray ? normalization_done : normalization_done[0];
}

function MatrixVectorProduct(matrixA, VectorB){
    // takes as input: a 2x2 squared matrix and a 2x1 column vector
    // returns: the product A.B as a p5.js vector object

    return createVector(matrixA[0][0]*VectorB.x + matrixA[0][1]*VectorB.y, 
                        matrixA[1][0]*VectorB.x + matrixA[1][1]*VectorB.y);
}

function MatrixTranspose(matrixA){
    // takes as input: a 2x2 squared matrix 
    // returns: the transposed matrix A^T

    return [[matrixA[0][0], matrixA[1][0]],
            [matrixA[0][1], matrixA[1][1]]];
}

function MatrixInverse(matrixA){
    // takes as input: a 2x2 squared matrix 
    // returns: the inverse matrix A^-1

    let det=matrixA[0][0]*matrixA[1][1] - matrixA[0][1]*matrixA[1][0];
    return [[matrixA[1][1], -matrixA[1][0]],
            [-matrixA[0][1], matrixA[0][0]]].map(row => row.map(val => val/det));
}

function MatrixMatrixProduct(matrixA, matrixB){
    // takes as input: two 2x2 squared matrix 
    // returns: the product A.B 

    return [[matrixA[0][0]*matrixB[0][0] + matrixA[0][1]*matrixB[1][0], matrixA[0][0]*matrixB[0][1] + matrixA[0][1]*matrixB[1][1]],
            [matrixA[1][0]*matrixB[0][0] + matrixA[1][1]*matrixB[1][0], matrixA[1][0]*matrixB[0][1] + matrixA[1][1]*matrixB[1][1]]];
}

function DoubleDotProduct(matrixA, matrixB){
    // takes as input: two 2x2 squared matrix 
    // returns: the double dot product A:B 
    let sum=0.0;
    for (let i = 0; i < 2; i++) {
      for (let j = 0; j < 2; j++) {
        sum += matrixA[i][j] * matrixB[i][j];
      }
    }
    return sum;
}

function CrossProduct(vecA, vecB){
  // takes as input: two 3D vectors
  // returns: the cross product 
  return createVector(vecA.y*vecB.z - vecA.z*vecB.y, vecA.z*vecB.x - vecA.x*vecB.z, vecA.x*vecB.y - vecA.y*vecB.x);
}

function HeavisideStepFunction(s, m=0.0){
  if (s<0.0) return 0.0;
  if (s>0.0) return 1.0;
  if (s==0.0) return m;
}

function LogisticFunction(t, param){
  // logistic function A/(1+exp(-w(t-tau))) with param=[A, w, tau]
  return param[0]/(1.0+exp(-param[1]*(t-param[2])));
}

function Rotation3DAroundAxis(N, u, alpha, to_normalize){
  // Rodrigues formula to rotate a 3D vector expressed in (x,y,z) around a unit-length axis N expressed also in (x,y,z) 
  // https://fr.wikipedia.org/wiki/Rotation_vectorielle

  u=u.copy().mult(cos(alpha)).add(N.copy().mult((1.0-cos(alpha))*u.dot(N))).add(CrossProduct(N, u).mult(sin(alpha)));
  
  if (to_normalize){
    return u.normalize();
  } else{
    return u;
  }
}




