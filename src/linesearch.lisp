(in-package #:matlisp-optimization)
(named-readtables:in-readtable :infix-dispatch-table)

(define-condition exceeded-maximum-iterations (warning)
  ((message :initarg :message))
  (:report (λ (c stream)
	      (when (slot-boundp c 'message)
		(format stream "~a~%" (slot-value c 'message))))))

(defun backtracking-linesearch (x0 g f df.g func &key (t0 1.0d0) (c 0.5) (rho 0.5) (max-iterations 20)
				&aux (x (zeros (dimensions x0) (class-of x0))) (nk 0))
  (tagbody start
     (iter (for k from 0 below max-iterations)
	   (for tk first t0 then (* tk rho))
	   (when (<= (funcall func (axpy! tk g (copy! x0 x)) 0) (+ f (* c tk df.g))) (finish))
	   (finally (setf (values t0 nk) (values tk (+ nk k)))
		    (when (= k max-iterations)
		      (restart-case (warn 'exceeded-maximum-iterations :message "Backtracking failed.")
			(continue-with-linesearch () (go start)))))))
  (values x t0 nk))

#+nil
(defun cauchy (delta g H)
  (let* ((H.g #i(H * g)) (gHg (dot g H.g))
	 (ng (norm g)) (plen (/ delta ng)))
    (scal! (if (<= gHg 0) plen (min (/ (* ng ng) gHg) plen)) (copy! g H.g))))
