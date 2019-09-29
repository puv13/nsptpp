# Customisation for Emacs CCmode 

Emacs CCmode is used to standardise the source code formatting throughout the project. 
We use Stroustrup's style with slight modifications. 

	;; CCmode customisation (add to .emacs)
	;; ----------------------------------------------------------------------
	

	;; New style defintion
	(c-add-style "nsptpp" 
		'("stroustrup"  ; based on 'stroustrup' style 
		(c-set-offset 'access-label '/) ; change indentation rules 
		(c-set-offset 'inclass '+)      ; indent public/privat keywords
	    (setq c-basic-offset 2)         ; default indent
		(setq fill-column 80)           ; 
		(auto-fill-mode 1)              ; break comment lines after 80 chars
		))


	;; Use this style for all languages
	(add-hook 'c-mode-common-hook '(lambda () 
		(c-set-style "nsptpp") 	      
		(setq comment-multi-line  t)    ; enable multi-line comments 
		(setq comment-continue  " * ")  ; start multi-line comments with *
		))	

	;; Language specific stuff for 'C'
	(add-hook 'c-mode-hook '(lambda () 
		(c-set-style "nsptpp")
		(setq comment-style 'multi-line)
		(setq comment-multi-line  t)
		(setq comment-continue  " * ")
		))	

	;; Language specific stuff for 'C++'
	(add-hook 'c++-mode-hook '(lambda () 
		(c-set-style "nsptpp") 	      
        (c-set-offset 'access-label '/) ; change indentation rules 
        (c-set-offset 'inclass '+)      ;
		(setq c-basic-offset 2)         ; default indent 

		))


