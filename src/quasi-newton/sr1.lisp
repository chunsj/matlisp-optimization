(in-package #:sr1)
(named-readtables:in-readtable :infix-dispatch-table)

;;Flip arguments for the inverse Hessian update.
(defun update! (Δx Δ∂f Bkp &optional (r 1d-6))
  ;;B^+ = B + ρ eta ⊗ eta ;;eta = Δ∂f - B Δx, ρ = 1/Δx·eta
  (let* ((eta (gem! -1 Bkp Δx 1 (copy Δ∂f)))
	 (eta.Δx (dot eta Δx)))
    (when (< (* r (norm eta) (norm Δx)) (abs eta.Δx))
      (ger! (/ eta.Δx) eta eta Bkp))
    Bkp))

;;Sparse SR1
(defun %lsr1-query (q0 buf &optional B0)
  (let ((q (etypecase B0
	     (null (copy q0))
	     (tensor-matrix (gem! 1 B0 q0 1 (copy q0)))
	     ((or number tensor-vector) (axpy! 1 q0 (scal B0 q0))))))
    (iter (for (eta . rho) in-dlist buf) (axpy! (* rho (dot eta q0)) eta q))
    q))

(defun %lsr1-update! (Δx Δ∂f buf &optional push? B0 (r 1d-6))
  (let* ((eta (axpy! 1 Δ∂f (scal! -1 (%lsr1-query Δx buf B0))))
	 (eta·Δx (dot eta Δx)))
    (when (< (* r (norm eta) (norm Δx)) (abs eta·Δx))
      (let* ((lbuf (if push? (dlist:dpush (cons nil nil) buf) (dlist:drdc buf)))
	     (ccon (dlist:dcar lbuf)))
	(setf (car ccon) eta (cdr ccon) (/ eta·Δx))
	lbuf))))

(defclass lsr1 ()
  ((buffer :initform nil) (B0 :initform nil) (eps :initarg :eps :initform 1d-6)
   (rank :initform 0) (n :initarg :n)))

(defun l-query (Δx lsr1)
  (declare (type lsr1 lsr1) (type (and dense-tensor tensor-vector) Δx))
  (%lsr1-query Δx #i(lsr1.buffer) #i(lsr1.B0)))

(defun l-update! (Δx Δ∂f lsr1)
  (declare (type lsr1 lsr1)
	   (type (and dense-tensor tensor-vector) Δx Δ∂f))
  (let ((buffer+ (%lsr1-update! Δx Δ∂f #i(lsr1.buffer) (if (< #i(lsr1.rank) #i(lsr1.n)) t) #i(lsr1.B0) #i(lsr1.eps))))
    (when buffer+
      (setf #i(lsr1.buffer) buffer+)
      (incf #i(lsr1.rank)))
    lsr1))
;;

#+nil
(let ((lsr1 (make-instance 'lsr1 :n 10))
      (A (randn '(10 10)))
      (x (randn 10)))
  (let ((y (randn 10))) (l-update! y #i(A * y) lsr1))
  (l-update! x #i(A * x) lsr1)
  (norm (t:- #i(A * x) (l-query x lsr1))))

#+nil
(let ((B (zeros '(10 10)))
      (A (randn '(10 10)))
      (x (randn 10)))
  (let ((y (randn 10))) (update! y #i(A * y) B))
  (update! x #i(A * x) B)
  (norm #i(A * x - B * x)))

#+nil
(let ((lsr1 (make-instance 'lsr1 :n 10))
      (A (randn '(10 10)))
      (x (randn 10)))
  ;;(setf #i(lsr1.buffer) (%lsr1-update! x #i(A * x) #i(lsr1.buffer) t))
  (lsr1-update! x #i(A * x) lsr1)
  (norm (tb- #i(A * x) (lsr1-query x lsr1)))
  #i (lsr1.rank)
  )

#+nil
(defmacro sr1-wrap (func)
  (binding-gensyms (gm gf)
    `(lambda (&rest args)
       ))
  )

#+nil
(defun sr1-descent (x0 func &key (atol 1d-6) (max-iterations 100))
  (let ((xk (copy x0)) (tk 1d0)
	(yk (zeros (dims x0))) (dk (zeros (dims x0)))
	(Hk (eye (make-list 2 :initial-element (dimensions x0 0))))
	(niters 0))
    (iter (for it from 0 below max-iterations)
	  (letv* ((f df (funcall func xk 1))))
	  (multiple-value-bind (f df) (funcall func xk 1)
	    (when (< (norm df) atol) (finish))
	    (when (> it 0) (sr1-update! tk dk (axpy! 1 df yk) Hk))
	    (gemv! -1 Hk df 0 dk)
	    (let ((g.d (dot dk df)))
	      (when (> g.d 0) (scal! -1 dk) (setf g.d (- g.d)))
	      (setq tk (nth-value 1 (backtracking-linesearch! xk dk f g.d func :t0 1d0 :c 0.1 :rho 0.5 :max-iterations 10))))
	    (copy! df yk) (scal! -1 yk))
	  (finally (incf niters it)
		   (when (= it max-iterations)
		     (warn 'exceeded-maximum-iterations :message "SR1 exceeded max-iterations."))))
    (values xk niters)))

#+nil
(defun sr1-descent (x0 func &key (atol 1d-6) (max-iterations 100))
  (let ((xk (copy x0)) (tk 1d0)
	(yk (zeros (dims x0))) (dk (zeros (dims x0)))
	(Sk (eye (make-list 2 :initial-element (dimensions x0 0))))
	(niters 0))
    (iter (for it from 0 below max-iterations)
	  (letv* ((f df (funcall func xk 1)))
	    (when (< (norm df) atol) (finish))
	    (when (> it 0) (sr1-update! dk (axpy! 1 df yk) Sk))
	    (gemv! -1 Hk df 0 dk)
	    
	    (let ((g.d (dot dk df)))
	      (when (> g.d 0) (scal! -1 dk) (setf g.d (- g.d)))
	      (setq tk (nth-value 1 (backtracking-linesearch! xk dk f g.d func :t0 1d0 :c 0.1 :rho 0.5 :max-iterations 10))))
	    (copy! df yk) (scal! -1 yk))
	  (finally (incf niters it)
		   (when (= it max-iterations)
		     (warn 'exceeded-maximum-iterations :message "SR1 exceeded max-iterations."))))
    (values xk niters)))
