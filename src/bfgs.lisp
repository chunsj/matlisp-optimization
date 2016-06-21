(in-package #:matlisp-optimization)
(named-readtables:in-readtable :infix-dispatch-table)

(defun rosenbrock (x &optional (grad 0))
  (m::compile-symbolic (double-float)
    (let* ((c (m::weylify '#i((1 - x[0])^2 + 100 * (x[1] - x[0]^2)^2)))
	   (dc (m::deriv c #स[x_0, x_1])))
      `(values-n (1+ grad) ,c ,dc ,(m::deriv dc #स[x_0, x_1])))))

(defun bfgs-descent (x.0 func
		     &rest keys
		     &key
		       (rtol 1d-6) (atol 1d-6) (max-iterations 100) (fd-δ 1d-3)
		     &allow-other-keys
		     &aux
		       (keytable (make-keytable keys))
		       (x.k (copy x.0)) (x.ls (zeros (dimensions x.0) (type-of x.0)))
		       (p.k (zeros (dimensions x.0) (type-of x.0)))
		       (δdf (zeros (dimensions x.0) (type-of x.0))) (δx (zeros (dimensions x.0) (type-of x.0)))
		       (H.k (eye (list (dimensions x.0 0) (dimensions x.0 0)))) (nk 0))
  (flet ((linefunc (tk &optional (grad 0))
	   (letv* ((f0 (funcall func (axpy! tk p.k (copy! x.k x.ls)) 0))
		   (fδ (if (< 0 grad) (funcall func (axpy! (* tk fd-δ) p.k x.ls) 0))))
	     (values-n (1+ grad) f0 (/ (- fδ f0) (* tk fd-δ))))))
    (tagbody start
       (iter (for it from 0 below max-iterations) (incf nk)
	     (letv* ((f.k df.k (funcall func x.k 1)))
	       (when (< (norm df.k) (+ atol (* rtol (norm x.k)))) (finish))
	       (when (< 0 nk) (bfgs:update! (axpy! 1 df.k δdf) (axpy! 1 x.k δx) H.k))
	       (gem! -1 H.k df.k 0 p.k)
	       (scal! -1 (copy! x.k δx)) (scal! -1 (copy! df.k δdf))
	       (letv* ((ti nk (bracketing-linesearch #'linefunc f.k (dot df.k p.k) keytable)))
		 (axpy! ti p.k x.k)))
	     (finally
	      (when (= it max-iterations)
		(restart-case (warn 'exceeded-maximum-iterations :message "LBFGS exceeded max-iterations.")
		  (continue-with-optimization? (answer) (when answer (go start)))))))))
  (values x.k nk))

(defun lbfgs-descent (x.0 func
		      &rest keys
		      &key
			(rtol 1d-6) (atol 1d-6) (max-iterations 100) (buffer-size 10) (fd-δ 1d-3)
		      &allow-other-keys
		      &aux
			(keytable (make-keytable keys))
			(x.k (copy x.0)) (x.ls (zeros (dimensions x.0) (type-of x.0)))
			(p.k (zeros (dimensions x.0) (type-of x.0)))
			(δdf (zeros (dimensions x.0) (type-of x.0))) (δx (zeros (dimensions x.0) (type-of x.0)))
			(H.k (make-instance 'bfgs:lbfgs :maximum-size buffer-size)) (nk 0))
  (flet ((linefunc (tk &optional (grad 0))
	   (letv* ((f0 (funcall func (axpy! tk p.k (copy! x.k x.ls)) 0))
		   (fδ (if (< 0 grad) (funcall func (axpy! (* tk fd-δ) p.k x.ls) 0))))
	     (values-n (1+ grad) f0 (/ (- fδ f0) (* tk fd-δ))))))
    (tagbody start
       (iter (for it from 0 below max-iterations) (incf nk)
	     (letv* ((f.k df.k (funcall func x.k 1)))
	       (when (< (norm df.k) (+ atol (* rtol (norm x.k)))) (finish))
	       (when (< 0 nk) (bfgs:l-update! (axpy! 1 df.k δdf) (axpy! 1 x.k δx) H.k))
	       (scal! -1 (copy! (bfgs:l-query df.k H.k) p.k))
	       (scal! -1 (copy! x.k δx)) (scal! -1 (copy! df.k δdf))
	       (letv* ((ti nk (bracketing-linesearch #'linefunc f.k (dot df.k p.k) keytable)))
		 (axpy! ti p.k x.k)))
	     (finally
	      (when (= it max-iterations)
		(restart-case (warn 'exceeded-maximum-iterations :message "LBFGS exceeded max-iterations.")
		  (continue-with-optimization? (answer) (when answer (go start)))))))))
  (values x.k nk))

#+nil

(defparameter *dbg* (randn 2))
;; (time (bfgs-descent *dbg* #'rosenbrock))
;; (time (lbfgs-descent *dbg* #'rosenbrock :buffer-size 10))
