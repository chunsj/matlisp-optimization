(in-package #:matlisp-optimization)
(named-readtables:in-readtable :infix-dispatch-table)

(define-condition exceeded-maximum-iterations (warning)
  ((message :initarg :message))
  (:report (λ (c stream)
	      (when (slot-boundp c 'message)
		(format stream "~a~%" (slot-value c 'message))))))
;;
(definline backtracking-update (tfn a-min a-max rho)
  (declare (ignore tfn))
  (+ a-min (* rho (- a-max a-min))))

;; [Nocedal & Wright, p.58]
(defun quadratic-update (tfn a-min a-max rho)
  (match tfn
    ((vector* (list (= 0) f0) (list (= 0) df0 1) (list t1 f1) _)
     (let ((a1 (/ (* -1 df0 t1 t1) (* 2 (- f1 f0 (* t1 df0))))))
       (if (<= a-min a1 a-max) a1
	   (backtracking-update tfn a-min a-max rho))))
    (_ (poly-update tfn a-min a-max rho 2))))

;; [Nocedal & Wright, p.58]
(defun cubic-update (tfn a-min a-max rho)
  (match tfn
    ((vector* (list (= 0) f0) (list (= 0) df0 1) (list t1 f1) (list t2 f2) _)
     (letv* ((δf1 (- f1 f0 (* t1 df0))) (δf2 (- f2 f0 (* t2 df0)))
	     (as (+ (* (expt t2 2) δf1) (* -1 (expt t1 2) δf2)))
	     (bs (+ (* -1 (expt t2 3) δf1) (* (expt t1 3) δf2)))
	     (det (* (expt (* t1 t2) 2) (- t1 t2)))
	     (bm (/ bs 3 as)) (dm (/ (sqrt (- (expt bs 2) (* 3 as det))) 3 (abs as)))
	     (a1 (- dm bm)) (a2 (- (- dm) bm)))
       (cond
	 ((and (realp a1) (<= a-min a1 a-max)) a1)
	 ((and (realp a2) (<= a-min a2 a-max)) a2)
	 (t (backtracking-update tfn a-min a-max rho)))))
    ((vector* (list (= 0) f0) (list (= 0) df0 1) (list t1 f1) (list (eql t1) df1 1) _)
     (letv* ((δdf (- df1 df0)) (δf (- f1 f0 (* t1 df0)))
	     (as (+ (* t1 δdf) (* -2 δf)))
	     (bs (+ (* -1 (expt t1 2) δdf) (* 3 t1 δf)))
	     (det (expt t1 3))
	     (bm (/ bs 3 as)) (dm (/ (sqrt (- (expt bs 2) (* 3 as det))) 3 (abs as)))
	     (a1 (- dm bm)) (a2 (- (- dm) bm)))
       (cond
	 ((and (realp a1) (<= a-min a1 a-max)) a1)
	 ((and (realp a2) (<= a-min a2 a-max)) a2)
	 (t (/ (+ a-min a-max) 2)))))
    (_ (poly-update tfn a-min a-max rho 3))))

(defun poly-update (tfn a-min a-max rho &optional n)
  (let* ((p (polyfit tfn n))
	 (as (roots (t:.* (range 1 (dimensions p 0) 1d0) (subtensor~ p '((1 nil)))))))
    (or (iter (for i from 0 below (dimensions as 0))
	      (match (ref as i)
		((complex (guard re (<= a-min re a-max)) (= 0)) (maximizing re))))
	(backtracking-update tfn a-min a-max rho))))
;;
(defmacro with-keytable (keys table &body body)
  (using-gensyms (decl (table) (vv ep))
    `(let* (,@decl
	    ,@(iter (for k in keys)
		    (letv* (((k &optional default) (ensure-list k))
			    (kk (intern (symbol-name k) :keyword)))
		      (collect
			  `(,k ,(if default
				    `(letv* ((,vv ,ep (if ,table (gethash ,kk ,table))))
				       (if ,ep ,vv ,default))
				    `(if ,table (gethash ,kk ,table))))))))
       ,@body)))

(defun make-keytable (key-values &aux (table (make-hash-table)))
  (iter (for (key val . _) on key-values by #'cddr)
	(setf (gethash key table) val))
  table)

(defun bracketing-linesearch (mfunc m.0 dm.0 &optional keytable)
  (with-keytable ((a-min 0d0) (a-max 1d0) (t.0 a-max) (c1 1d-4) (c2 0.5d0) (rho 0.5d0) (linesearch-max-iterations 20)) keytable
    (let* ((nk 0) (tfn (make-array 2 :initial-contents `((0d0 ,m.0) (0d0 ,dm.0 1)) :adjustable t :fill-pointer t))
	   (m.a-min (if (= a-min 0) m.0 (funcall mfunc a-min 0))) (m.a-max (funcall mfunc a-max 0)))
      (macrolet ((optimization-loop (&body kernel)
		   `(tagbody start
		       (iter (for k from 0 below linesearch-max-iterations)
			     (for t.k first (min a-max t.0) then (case (length tfn)
								   (2 (opt:backtracking-update tfn a-min a-max rho))
								   (3 (opt:quadratic-update tfn a-min a-max rho))
								   (t (opt:cubic-update tfn a-min a-max rho))))
			     (progn ,@kernel)
			     (finally (setf (values t.0 nk) (values t.k (+ nk k 1)))
				      (when (= k linesearch-max-iterations)
					(restart-case (warn 'opt:exceeded-maximum-iterations :message "Backtracking failed.")
					  (continue-with-linesearch () (go start)))))))))
	(if c2
	    (optimization-loop
	     (letv* ((m.k dm.k error? (handler-case (funcall mfunc t.k 1)
					(t () (values m.0 nil t)))))
	       (setf (fill-pointer tfn) 2)
	       (if error? (continue)
		   (progn
		     (vector-push-extend (list t.k m.k) tfn)
		     (vector-push-extend (list t.k dm.k 1) tfn)))
	       ;;Feasible point
	       (when (and (<= m.k (+ m.0 (* c1 t.k dm.0))) ;;Wolfe-1
			  (<= (* c2 dm.0) dm.k 0))         ;;Wolfe-2
		 (finish))
	       ;;update bounds [Mor\'e Thuente '94 ?]
	       (if (or (< m.a-min m.k) (< 0 dm.k))
		   (setf a-max t.k m.a-max m.k)
		   (setf a-min t.k m.a-min m.k))))
	    (optimization-loop
	     (letv* ((m.k dm.k error? (handler-case (funcall mfunc t.k 0)
					(t () (values m.0 nil t)))
			  :type real nil t))
	       (when (< 3 (length tfn))
		 (setf (aref tfn 2) (aref tfn (1- (length tfn))))
		 (setf (fill-pointer tfn) 3))
	       (if error? (continue)
		   (vector-push-extend (list t.k m.k) tfn))
	       ;;Feasible point
	       (when (<= m.k (+ m.0 (* c1 t.k dm.0)))     ;;Wolfe-1
		 (finish))
	       ;;update bounds
	       (if (< m.a-min m.k)
		   (setf a-max t.k m.a-max m.k))))))
      (values t.0 nk))))

;;
;; (range 1 10)
#+nil
(defun cauchy (delta g H)
  (let* ((H.g #i(H * g)) (gHg (dot g H.g))
	 (ng (norm g)) (plen (/ delta ng)))
    (scal! (if (<= gHg 0) plen (min (/ (* ng ng) gHg) plen)) (copy! g H.g))))
