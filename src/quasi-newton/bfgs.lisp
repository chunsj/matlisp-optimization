(in-package #:bfgs)
(named-readtables:in-readtable :infix-dispatch-table)

(define-condition bfgs-warning (warning)
  ((message :initarg :message))
  (:report (λ (c stream) (when (slot-boundp c 'message)	(format stream "~a~%" (slot-value c 'message))))))

;;Flip Δx Δ∂f for inverse Hessian update. The default is infact the DFP update.
(defun update! (Δx Δ∂f Bkp &optional (r 1d-6) &aux (Δx.Δ∂f (dot Δx Δ∂f)))
  ;;B^+ = P' * B * P + rho .* Δ∂f ⊗ Δ∂f
  ;;rho = 1/(Δx · Δ∂f)
  ;;P = (id - rho .* Δx ⊗ Δ∂f)
  (tagbody
   main
     (restart-case (unless (< r Δx.Δ∂f) (warn 'bfgs-warning :message "non positive-definite update. skipping...") (go end))
       (continue () (go bfgs)))
   bfgs
     (let ((rho (/ Δx.Δ∂f)) (B·Δx (gem 1 Bkp Δx nil nil)))
       (ger! (* (1+ (* (dot Δx B·Δx) rho)) rho) Δ∂f Δ∂f (ger! (- rho) Δ∂f B·Δx (ger! (- rho) B·Δx Δ∂f Bkp))))
   end)
  Bkp)

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
  (tagbody
   main
     (restart-case (unless (< r Δx·Δ∂f) (warn 'bfgs-warning :message "non positive-definite update. skipping...") (go end))
       (continue () (go bfgs)))
   bfgs
     (let* ((lbuf (if push? (dlist:dpush (list (zeros (dimensions Δx) (class-of Δx)) (zeros (dimensions Δ∂f) (class-of Δ∂f))) buf) (dlist:drdc buf)))
	    (bcon (dlist:dcar lbuf)))
       (copy! Δx (first bcon)) (copy! Δ∂f (second bcon))
       (setf (cddr bcon) (/ Δx·Δ∂f)
	     buf lbuf))
   end)
  buf)

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

;;Tests
(5am:test bfgs-test
  (let ((B (eye '(10 10)))
	(A (psd-proj (randn '(10 10))))
	(x (randn 10)))
    (handler-bind ((bfgs-warning #'(lambda (c) (invoke-restart 'continue))))
      (let ((y (randn 10))) (update! y #i(A * y) B))
      (update! x #i(A * x) B)
      (5am:is (< (norm #i(A * x - B * x)) (* 100 double-float-epsilon))))))

(5am:test lbfgs-test
  (let ((lbfgs (make-instance 'lbfgs :n 10))
	(A (psd-proj (randn '(10 10))))
	(x (randn 10)))
    (handler-bind ((bfgs-warning #'(lambda (c) (invoke-restart 'continue))))
      (let ((y (randn 10))) (l-update! y #i(A * y) lbfgs))
      (l-update! x #i(A * x) lbfgs)
      (5am:is (< (norm (t:- #i(A * x) (l-query x lbfgs))) (* 100 double-float-epsilon))))))
