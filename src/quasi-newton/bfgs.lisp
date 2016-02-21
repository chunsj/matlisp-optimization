(in-package #:bfgs)
(named-readtables:in-readtable :infix-dispatch-table)

(define-condition bfgs-error (error)
  ((message :initarg :message))
  (:report (λ (c stream) (when (slot-boundp c 'message)	(format stream "~a~%" (slot-value c 'message))))))

;;Flip Δx Δ∂f for Hessian update.
(defun update! (Δx Δ∂f Bkp &optional (r 1d-6) &aux (Δx.Δ∂f (dot Δx Δ∂f)))
  ;;B^+ = P' * B * P + rho .* Δ∂f ⊗ Δ∂f
  ;;rho = 1/(Δx · Δ∂f)
  ;;P = (id - rho .* Δx ⊗ Δ∂f)
  (assert (< r Δx.Δ∂f) nil 'bfgs-error)
  (let ((rho (/ Δx.Δ∂f)) (B·Δx (gem 1 Bkp Δx nil nil)))
    (ger! (* (1+ (* (dot Δx B·Δx) rho)) rho) Δ∂f Δ∂f (ger! (- rho) Δ∂f B·Δx (ger! (- rho) B·Δx Δ∂f Bkp)))))

;;
(defun %lbfgs-query (q0 buf &optional B0 &aux (q (copy q0)))
  (let ((alpha nil))
    (iter (for (y.k s.k . rho.k) in-dlist buf)
	  (let ((a.k (* rho.k (dot s.k q))))
	    (push a.k alpha)
	    (axpy! (- a.k) y.k q)))
    (let ((r (etypecase B0
	       (null q)
	       (tensor-matrix (gem 1 B0 q nil nil))
	       ((or number tensor-vector) (axpy! 1 q (scal! B0 q))))))
      (iter (for a.k in alpha)
	    (for (y.k s.k . rho.k) in-dlist (first buf) in-reverse t)
	    (axpy! (- a.k (* rho.k (dot y.k r))) s.k r))
      q)))

(defun %lbfgs-update! (Δx Δ∂f buf &optional push? (r 1d-6) &aux (Δx·Δ∂f (dot Δx Δ∂f)))
  (assert (< r Δx·Δ∂f) nil 'bfgs-error)
  (let* ((lbuf (if push? (dlist:dpush (list (zeros (dimensions Δx) (class-of Δx)) (zeros (dimensions Δ∂f) (class-of Δ∂f))) buf) (dlist:drdc buf)))
	 (bcon (dlist:dcar lbuf)))
    (copy! Δx (first bcon)) (copy! Δ∂f (second bcon))
    (setf (cddr bcon) (/ Δx·Δ∂f))
    lbuf))

(defclass lbfgs ()
  ((buffer :initform nil) (B0 :initform nil) (eps :initarg :eps :initform 1d-6)
   (rank :initform 0) (n :initarg :n)))

(defun l-query (Δx lbfgs)
  (declare (type lbfgs lbfgs) (type (and dense-tensor tensor-vector) Δx))
  (%lbfgs-query Δx #i(lbfgs.buffer) #i(lbfgs.B0)))

(defun l-update! (Δx Δ∂f lbfgs)
  (declare (type lbfgs lbfgs)
	   (type (and dense-tensor tensor-vector) Δx Δ∂f))
  (let ((buffer+ (%lbfgs-update! Δx Δ∂f #i(lbfgs.buffer) (if (< #i(lbfgs.rank) #i(lbfgs.n)) t) #i(lbfgs.eps))))
    (when buffer+
      (setf #i(lbfgs.buffer) buffer+)
      (incf #i(lbfgs.rank)))
    lbfgs))

#+nil
(let ((lbfgs (make-instance 'lbfgs :n 10))
      (A (randn '(10 10)))
      (x (randn 10)))
  (let ((y (randn 10))) (l-update! y #i(A * y) lbfgs))
  (l-update! x #i(A * x) lbfgs)
  (norm (t:- #i(A * x) (l-query x lbfgs))))

#+nil
(let ((B (eye '(10 10)))
      (A (psd-proj (randn '(10 10))))
      (x (randn 10)))
  (let ((y (randn 10))) (update! y #i(A * y) B))
  (update! x #i(A * x) B)
  (norm #i(A * x - B * x)))

#+nil
(defun lbfgs-descent (x0 func &key (atol 1d-6) (max-iterations 100) (buf-size 10))
  (let ((xk (copy x0)) (tk 1d0)
	(yk (zeros (dims x0))) (dk (zeros (dims x0)))
	(buf nil)
	(buf-count 0))
    (values xk
	    (iter (for it from 0 below max-iterations)
		  (multiple-value-bind (f df) (funcall func xk 1)
		    (print (list f (norm df) tk))
		    (when (< (norm df) atol) (finish))
		    (when (> it 0)
		      (setf buf (lbfgs-update! (scal! tk dk) (axpy! 1 df yk) buf (when (< buf-count buf-size) (incf buf-count)))))
		    (lbfgs-query! (progn (copy! df dk) (scal! -1 dk)) buf)
		    (setq tk (nth-value 1 (backtracking-linesearch! xk dk f (dot dk df) func :t0 1d0 :c 0.1 :rho 0.5 :max-iterations 25)))
		    (copy! df yk) (scal! -1 yk))
		  (finally (when (= it max-iterations)
			     (format t "Exceeded max-iterations.~%"))
			   (return it))))))

#+nil
(with-optimization (:debug 3)
  (defun bfgs-descent (x0 func &key (atol 1d-6) (max-iterations 100))
    (let ((xk (copy x0)) (tk 1d0)
	  (yk (zeros (dims x0))) (dk (zeros (dims x0)))
	  (Hk (eye (make-list 2 :initial-element (aref (dimensions x0) 0))))
	  (niters 0))
      (tagbody start
	 (iter (for it from 0 below max-iterations)
	       (multiple-value-bind (f df) (funcall func xk 1)
		 ;;(print (list f (norm df)))
		 (when (< (norm df) atol) (finish))
		 (when (> it 0) (bfgs-update! tk dk (axpy! 1 df yk) Hk))
		 (setq tk (nth-value 1 (backtracking-linesearch! xk (gemv! -1 Hk df 0 dk) f (dot dk df) func :t0 1d0 :c 0.1 :rho 0.5 :max-iterations 10)))
		 (copy! df yk) (scal! -1 yk))
	       (finally	(incf niters it)
			(when (= it max-iterations)
			  (restart-case (warn 'exceeded-maximum-iterations :message "BFGS exceeded max-iterations.")
			    (continue-with-optimization? (answer) (when answer (go start))))))))
      (values xk niters))))
