{{ fullname | escape | underline }}

.. automodule:: {{ fullname }}

    {% block classes %}
    {% if classes %}
    {% for item in classes %}
        .. autoclass:: {{ item }}
            :members:
            :undoc-members:
            :show-inheritance:
            :special-members: __init__

            .. autoclasstoc::
    {%- endfor %}
    {% endif %}
    {% endblock %}

    {% block functions %}
    {% if functions %}
    {% for item in functions %}
        .. autofunction:: {{ item }}
    {%- endfor %}
    {% endif %}
    {% endblock %}

    {% block attributes %}
    {% if attributes %}
    {% for item in attributes %}
        .. autoattribute:: {{ item }}
    {%- endfor %}
    {% endif %}
    {% endblock %}

    {% block exceptions %}
    {% if exceptions %}
    {% for item in exceptions %}
        .. autoexception:: {{ item }}
    {%- endfor %}
    {% endif %}
    {% endblock %}