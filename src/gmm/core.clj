(ns gmm.core
  (:use [incanter.core]))

(defn gaussian
  ; computes 2-dimensional gaussian pdf (as function of point (2-vector))
  ; μ is a 2-vector
  ; Σ is a 2x2 matrix
  [μ Σ]
  (fn [x]
    (*
      (Math/pow (* 2 (. Math PI)) -1)
      (det Σ)
      (Math/pow
        (. Math E)
        (mult
          -0.5
          (mmult
            (minus x μ)
            (solve Σ)
            (trans (minus x μ))))))))

(defn expectation
  ; returns responsibility matrix
  [points means covariances coeffs]
  (let
    [distributions (map gaussian means covariances)]
    (matrix
      (map
        (fn [dist coeff]
          (map
            (fn [point]
              (/
                (* coeff (dist point))
                (sum
                  (map
                    (fn [dist coeff] (* coeff (dist point)))
                    distributions
                    coeffs))))
            points))
        distributions
        coeffs))))

(defn maximization
  ; returns means, covariances, coeffs
  ; γ is the responsibility matrix
  [points γ]
  (let
    [N      (map sum γ)
     means  (map
              (fn [γ_k N_k]
                (div
                  (vec
                    (reduce plus
                      (map mult (map first points) γ_k)))
                  N_k))
              γ
              N)
     covars (map
              (fn [γ_k μ_k N_k]
                (div
                  (reduce plus
                    (map
                      (fn
                        [x_i y_ik]
                        (mult
                          y_ik
                          (mmult
                            (minus x_i μ_k)
                            (trans (minus x_i μ_k)))))
                      (map first points)
                      γ_k))
                  N_k))
              γ
              (map first means)
              N)
     coeffs (map
              #(/ % (nrow γ))
              N)
     ]
    [(map #(vector (vec %)) means) ; ((. .) (. .)) -> [[[. .]] [[. .]]]
     covars
     coeffs]))

(def points [[[0 0]] [[1 1]] [[2 2]]])
(def means [[[2 5]] [[-1 -2]]])
(def covariances [(matrix [[1 0] [0 1]]) (matrix [[1 0] [0 1]])])
(def coeffs [1 1])

; (defn expm
;   ; computes matrix exponential using JBLAS
;   ; (converts to and from JBLAS representation)
;   [oldMat]
;   (let
;     [newMat
;       (DoubleMatrix.
;         (java.util.ArrayList.
;           (flatten
;             (to-list oldMat))))]
;     (.reshape newMat 2 2)
;     (matrix
;       (.toArray
;        (MatrixFunctions/expm newMat))
;       3)))
