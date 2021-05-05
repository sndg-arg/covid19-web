class JobValidationError(Exception):
    def __init__(self, status_text,errors=None):
        self.status_text = status_text

        self.errors = errors if errors else []
        super().__init__(self.status_text)
