(in-package #:sr1)
(named-readtables:in-readtable :infix-dispatch-table)

(define-condition sr1-warning (warning)
  ((message :initarg :message))
  (:report (λ (c stream) (when (slot-boundp c 'message)	(format stream "~a~%" (slot-value c 'message))))))

;;Flip arguments for the inverse Hessian update.
(defun update! (Δx Δ∂f Bkp &optional (r 1d-6))
  ;;B^+ = B + ρ eta ⊗ eta ;;eta = Δ∂f - B Δx, ρ = 1/Δx·eta
  (let* ((eta (gem! -1 Bkp Δx 1 (copy Δ∂f)))
	 (eta.Δx (dot eta Δx)))
    (tagbody main
       (restart-case (unless (< (* r (norm eta) (norm Δx)) (abs eta.Δx)) (warn 'sr1-warning :message "orthogonal update. skipping...") (go end))
	 (continue () (go sr1)))
     sr1 (ger! (/ eta.Δx) eta eta Bkp)
     end)
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
    (tagbody main
       (restart-case (unless (< (* r (norm eta) (norm Δx)) (abs eta·Δx)) (warn 'sr1-warning :message "orthogonal update. skipping...") (go end))
	 (continue () (go sr1)))
     sr1
       (let* ((lbuf (if push? (dlist:dpush (cons nil nil) buf) (dlist:drdc buf)))
	      (ccon (dlist:dcar lbuf)))
	 (setf (car ccon) eta (cdr ccon) (/ eta·Δx)
	       buf lbuf))
     end)
    buf))

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

;;Tests
(5am:test sr1-test
  (let ((B (eye '(10 10)))
	(A (psd-proj (randn '(10 10))))
	(x (randn 10)))
    (handler-bind ((sr1-warning #'(lambda (c) (invoke-restart 'continue))))
      (let ((y (randn 10))) (update! y #i(A * y) B))
      (update! x #i(A * x) B)
      (5am:is (< (norm #i(A * x - B * x)) (* 100 double-float-epsilon))))))

(5am:test lsr1-test
  (let ((lsr1 (make-instance 'lsr1 :n 10))
	(A (psd-proj (randn '(10 10))))
	(x (randn 10)))
    (handler-bind ((sr1-warning #'(lambda (c) (invoke-restart 'continue))))
      (let ((y (randn 10))) (l-update! y #i(A * y) lsr1))
      (l-update! x #i(A * x) lsr1)
      (5am:is (< (norm (t:- #i(A * x) (l-query x lsr1))) (* 100 double-float-epsilon))))))
