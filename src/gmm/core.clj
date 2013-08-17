(ns gmm.core
  (:use [incanter core stats]))

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

(defn gmm
  ([points n]
    (let
      [xs       (map #(first (first %)) points)
       ys       (map #(second (first %)) points)
       x-coords (repeatedly n #(+ (apply min xs)
                                  (* (rand) (- (apply max xs) (apply min xs)))))
       y-coords (repeatedly n #(+ (apply min ys)
                                  (* (rand) (- (apply max ys) (apply min ys)))))
       means    (map #(vector (vector %1 %2)) x-coords y-coords)
       covars   (repeat n (matrix [[1 0] [0 1]]))
       coeffs   (repeat n 1)]
      (gmm points means covars coeffs)))
  ([points means covars coeffs]
    (println means)
    (let
      [γ (expectation points means covars coeffs)
       [new-means new-covars new-coeffs] (maximization points γ)]
      (gmm points new-means new-covars new-coeffs))))
