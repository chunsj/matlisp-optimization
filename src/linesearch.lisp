(in-package #:matlisp-optimization)
(named-readtables:in-readtable :infix-dispatch-table)

(define-condition exceeded-maximum-iterations (warning)
  ((message :initarg :message))
  (:report (λ (c stream)
	      (when (slot-boundp c 'message)
		(format stream "~a~%" (slot-value c 'message))))))

(defun backtracking-update (tfn rho)
  (ematch (aref tfn (1- (length tfn)))
    ((λlist tk _ &optional _) (* tk rho))))

;; [Nocedal & Wright, p.58]
(defun quadratic-update (tfn)
  (match tfn
    ((vector* (list (= 0) f0) (list (= 0) df0 1) (list t1 f1) _)
     (/ (* -1 df0 t1 t1) (* 2 (- f1 f0 (* t1 df0)))))
    (_ (poly-update tfn 0 1 2))))

;; [Nocedal & Wright, p.58]
(defun cubic-update (tfn a-min a-max)
  (match tfn
    ((vector* (list (= 0) f0) (list (= 0) df0 1) (list t1 f1) (list t2 f2) _)
     (letv* ((δf1 (- f1 f0 (* t1 df0))) (δf2 (- f2 f0 (* t2 df0)))
	     (as (+ (* (expt t2 2) δf1) (* -1 (expt t1 2) δf2)))
	     (bs (+ (* -1 (expt t2 3) δf1) (* (expt t1 3) δf2)))
	     (det (* (expt (* t1 t2) 2) (- t1 t2)))
	     (bm (/ bs 3 as)) (dm (/ (sqrt (- (expt bs 2) (* 3 as det))) 3 (abs as)))
	     (a1 (- dm bm)) (a2 (- (- dm) bm)))
       (cond
	 ((<= a-min a1 a-max) a1)
	 ((<= a-min a2 a-max) a2)
	 (t (/ (+ a-min a-max) 2)))))
    ((vector* (list (= 0) f0) (list (= 0) df0 1) (list t1 f1) (list (eql t1) df1 1) _)
     (letv* ((δdf (- df1 df0)) (δf (- f1 f0 (* t1 df0)))
	     (as (+ (* t1 δdf) (* -2 δf)))
	     (bs (+ (* -1 (expt t1 2) δdf) (* 3 t1 δf)))
	     (det (expt t1 3))
	     (bm (/ bs 3 as)) (dm (/ (sqrt (- (expt bs 2) (* 3 as det))) 3 (abs as)))
	     (a1 (- dm bm)) (a2 (- (- dm) bm)))
       (cond
	 ((<= a-min a1 a-max) a1)
	 ((<= a-min a2 a-max) a2)
	 (t (/ (+ a-min a-max) 2)))))
    (_ (poly-update tfn a-min a-max 3))))

(defun poly-update (tfn a-min a-max &optional n)
  (let* ((p (polyfit tfn n))
	 (as (roots (t:.* (range 1 (dimensions p 0) 1d0) (subtensor~ p '((1 nil)))))))
    (or (iter (for i from 0 below (dimensions as 0))
	      (match (ref as i)
		((complex (guard re (<= a-min re a-max)) (= 0)) (maximizing re))))
	(/ (+ a-min a-max) 2))))
;;
(defun backtracking-linesearch (f f0 df0 &key (t0 1.0d0) (c 0.5) (rho 0.5) (max-iterations 20) &aux (nk 0))
  (declare (type function f))
  (tagbody start
     (iter (for k from 0 below max-iterations)
	   (for tk first t0 then (* tk rho))
	   (when (<= (funcall f tk) (+ f0 (* c tk df0))) (finish))
	   (finally (setf (values t0 nk) (values tk (+ nk k)))
	     (when (= k max-iterations)
	       (restart-case (warn 'exceeded-maximum-iterations :message "Backtracking failed.")
		 (continue-with-linesearch () (go start)))))))
  (values t0 nk))

#+nil
(let ((fval #((0d0 0.1d0) (0d0 -0.1d0 1) (0.5d0 0.2d0) (0.5d0 0.3d0 1))))
  (list (quadratic-update fval)
	(poly-update fval 0.0 0.7 2))
  )


;; (range 1 10)

#+nil
(defun cauchy (delta g H)
  (let* ((H.g #i(H * g)) (gHg (dot g H.g))
	 (ng (norm g)) (plen (/ delta ng)))
    (scal! (if (<= gHg 0) plen (min (/ (* ng ng) gHg) plen)) (copy! g H.g))))
